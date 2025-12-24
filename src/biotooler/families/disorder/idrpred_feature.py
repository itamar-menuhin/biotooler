"""IDRPred backend for disorder prediction.

This module provides protein intrinsic disorder prediction using IDRPred.
"""

import subprocess
from io import StringIO

import numpy as np
from Bio.SeqRecord import SeqRecord

from biotooler.core.record import get_molecule_type
from biotooler.core.translation import ensure_protein_record
from biotooler.families.disorder.integration import require_idrpred_cli
from biotooler.features.aggregation import AggregationSpec, PositionSpace


class IDRPredConsensusMask:
    """Compute per-residue IDR membership mask using IDRPred consensus predictor.

    This feature implements the PositionalFeature protocol to compute protein
    intrinsic disorder regions (IDRs) using IDRPred's consensus prediction method.
    It returns a binary mask indicating which residues are predicted to be in
    disordered regions.

    The feature operates in RESIDUE position space (per amino acid) and returns a
    binary vector (IDRPRED_IDR values) for each position in the protein sequence,
    where 1 indicates the residue is in a predicted IDR and 0 indicates it is not.

    Input handling:
    - Protein sequences: Used directly for prediction
    - DNA/RNA sequences: Automatically translated to protein using ensure_protein_record()

    Settings:
    - Uses IDRPred CLI tool via subprocess
    - Returns binary mask (0 or 1) for IDR membership
    - Deterministic: same input + same config → same mask

    Args:
        table: NCBI genetic code table number for translation (default: 1)
        use_orf_if_present: If True, use attached ORF for translation (default: True)
        strip_terminal_stop: If True, strip terminal stop codon (default: True)
        on_internal_stop: Action for internal stops: "error" (raise exception) or
            "ignore" (continue translation) (default: "error")

    Examples:
        >>> from Bio.Seq import Seq
        >>> from Bio.SeqRecord import SeqRecord
        >>> feature = IDRPredConsensusMask()
        >>> # Protein sequence
        >>> record = SeqRecord(Seq("MKALVSWGR"), id="test")
        >>> record.annotations["molecule_type"] = "protein"
        >>> result = feature.compute_vector(record)
        >>> print(result["IDRPRED_IDR"])  # Binary array of IDR membership
        >>> # DNA sequence (will be translated)
        >>> dna = SeqRecord(Seq("ATGAAAGCCCTGGTG"), id="test")
        >>> dna.annotations["molecule_type"] = "DNA"
        >>> result = feature.compute_vector(dna)
        >>> len(result["IDRPRED_IDR"])  # Length matches translated protein

    Notes:
        - idrpred CLI is checked at feature usage time
        - Results are deterministic for the same input sequence
        - Empty sequences return empty arrays
        - Translation follows the same rules as other protein family features
        - Output key IDRPRED_IDR is distinct from metapredict's DISORDER_P
    """

    def __init__(
        self,
        *,
        table: int = 1,
        use_orf_if_present: bool = True,
        strip_terminal_stop: bool = True,
        on_internal_stop: str = "error",
    ):
        """Initialize IDRPredConsensusMask.

        Args:
            table: NCBI genetic code table for translation (default: 1)
            use_orf_if_present: Use attached ORF if present (default: True)
            strip_terminal_stop: Strip terminal stop codon (default: True)
            on_internal_stop: "error" (raise exception) or "ignore" (continue
                translation) for internal stops (default: "error")
        """
        self.table = table
        self.use_orf_if_present = use_orf_if_present
        self.strip_terminal_stop = strip_terminal_stop
        self.on_internal_stop = on_internal_stop

    @property
    def position_space(self) -> PositionSpace:
        """Return RESIDUE position space.

        IDR prediction is computed per amino acid residue.

        Returns:
            PositionSpace.RESIDUE
        """
        return PositionSpace.RESIDUE

    @property
    def vector_keys(self) -> dict[str, AggregationSpec]:
        """Return aggregation specifications for each feature key.

        The IDRPRED_IDR (IDR membership) values are aggregated using mean,
        which gives the fraction of residues in IDRs over a window.

        Returns:
            Dictionary mapping "IDRPRED_IDR" to AggregationSpec with np.mean aggregation
        """
        return {
            "IDRPRED_IDR": AggregationSpec(name="MEAN", aggregation_fn=np.mean),
        }

    def compute_vector(self, record: SeqRecord, **kwargs) -> dict[str, np.ndarray]:
        """Compute per-residue IDR membership mask for the entire sequence.

        This method uses IDRPred to predict intrinsic disorder regions
        for each amino acid position. If the input is DNA/RNA, it is automatically
        translated to protein using ensure_protein_record().

        The computation uses IDRPred's CLI tool with FASTA input/output via stdin/stdout.

        Args:
            record: SeqRecord containing the DNA, RNA, or protein sequence to analyze
            **kwargs: Additional parameters that may include:
                - positions: Optional array of 0-based position indices. If provided,
                  the feature still computes the full vector (for correct context),
                  and the returned array will have the full length with valid values
                  at all positions. The window engine will extract needed positions
                  from the full vector.

        Returns:
            Dictionary with key "IDRPRED_IDR" mapping to numpy array of binary IDR
            membership values. Array length equals the protein sequence length, with
            values in {0, 1}.

        Raises:
            ImportError: If idrpred CLI is not available on PATH
            ValueError: If translation fails or other errors occur
            subprocess.CalledProcessError: If idrpred command fails

        Notes:
            - Uses IDRPred's consensus prediction method
            - Automatically handles DNA/RNA to protein translation
            - Returns deterministic results for the same input sequence
            - Computation is performed on the full sequence for correct context
            - positions parameter is accepted but full computation is always performed
        """
        # Check for idrpred CLI availability
        require_idrpred_cli()

        # Get molecule type to check if translation is needed
        mol_type = get_molecule_type(record)

        # If DNA or RNA, translate to protein
        if mol_type.upper() in ("DNA", "RNA"):
            protein_record = ensure_protein_record(
                record,
                table=self.table,
                use_orf_if_present=self.use_orf_if_present,
                strip_terminal_stop=self.strip_terminal_stop,
                on_internal_stop=self.on_internal_stop,
            )
        else:
            # Already protein, use as-is
            protein_record = record

        # Get protein sequence string
        protein_seq = str(protein_record.seq)

        # Handle empty sequence
        if len(protein_seq) == 0:
            return {"IDRPRED_IDR": np.array([], dtype=np.float64)}

        # Prepare FASTA input for idrpred
        fasta_input = f">seq\n{protein_seq}\n"

        # Run idrpred via subprocess with stdin/stdout
        try:
            result = subprocess.run(
                ["idrpred", "-", "-"],
                input=fasta_input,
                capture_output=True,
                text=True,
                check=True,
            )
            tsv_output = result.stdout
        except subprocess.CalledProcessError as e:
            raise ValueError(
                f"IDRPred command failed for record {record.id!r}: {e.stderr}"
            ) from e
        except FileNotFoundError as e:
            # This should not happen due to require_idrpred_cli, but handle it anyway
            raise ImportError(
                "idrpred command not found. Ensure it is installed and on PATH."
            ) from e

        # Parse TSV output to get IDR regions
        idr_regions = parse_tsv_regions(tsv_output)

        # Convert regions to per-residue mask
        idr_mask = regions_to_mask(idr_regions, len(protein_seq))

        return {"IDRPRED_IDR": idr_mask}


def parse_tsv_regions(tsv_output: str) -> list[tuple[int, int]]:
    """Parse IDRPred TSV output to extract IDR regions.

    This parser is robust to:
    - Empty lines
    - Comment/header lines (lines starting with #)
    - Extra columns beyond start/end coordinates
    - Different coordinate systems (defaults to 1-based inclusive)

    Args:
        tsv_output: TSV output from idrpred command

    Returns:
        List of (start, end) tuples representing IDR regions in 1-based inclusive
        coordinates. Empty list if no IDRs are predicted.

    Notes:
        - Assumes 1-based inclusive coordinates by default
        - If header explicitly indicates 0-based, will handle that (future)
        - Ignores empty lines and comment lines
        - Extracts first two integer columns as start/end
        - Deterministic: same input → same output
    """
    regions = []

    # Read TSV line by line
    for line in StringIO(tsv_output):
        line = line.strip()

        # Skip empty lines
        if not line:
            continue

        # Skip comment/header lines
        if line.startswith("#"):
            continue

        # Split by tab and extract first two columns as integers
        parts = line.split("\t")
        if len(parts) < 2:
            continue

        try:
            # Extract first two columns as start and end
            # IDRPred outputs 1-based inclusive coordinates
            start = int(parts[0])
            end = int(parts[1])

            # Validate and add region
            if start > 0 and end > 0 and start <= end:
                regions.append((start, end))
        except (ValueError, IndexError):
            # Skip lines that don't have valid integer coordinates
            continue

    return regions


def regions_to_mask(regions: list[tuple[int, int]], seq_length: int) -> np.ndarray:
    """Convert IDR regions to a per-residue binary mask.

    Args:
        regions: List of (start, end) tuples in 1-based inclusive coordinates
        seq_length: Length of the protein sequence

    Returns:
        Binary numpy array of shape (seq_length,) where 1 indicates the residue
        is in an IDR and 0 indicates it is not. dtype is float64 for consistency
        with other features.

    Notes:
        - Coordinates are clamped to [1, seq_length]
        - 1-based inclusive coordinates are converted to 0-based indices
        - Overlapping regions are merged (each position is 0 or 1)
        - Empty regions list returns all zeros
    """
    # Initialize mask with zeros
    mask = np.zeros(seq_length, dtype=np.float64)

    # Fill in IDR regions
    for start, end in regions:
        # Clamp coordinates to valid range [1, seq_length]
        start = max(1, min(start, seq_length))
        end = max(1, min(end, seq_length))

        # Skip invalid regions
        if start > end:
            continue

        # Convert to 0-based indices and mark region as IDR
        # 1-based inclusive [start, end] -> 0-based [start-1, end)
        mask[start - 1 : end] = 1.0

    return mask

"""Reference sequence set for storing and processing CDS and protein sequences."""

from pathlib import Path
from typing import Self

from Bio import SeqIO
from Bio.Seq import Seq


class ReferenceSequenceSet:
    """Container for reference CDS and protein sequences with normalization and translation.

    Stores CDS sequences (DNA/RNA) and optionally protein sequences.
    Provides normalized access to sequences with validation and caching.

    Attributes:
        cds: Dictionary mapping sequence IDs to CDS sequences (DNA/RNA)
        proteins: Optional dictionary mapping sequence IDs to protein sequences
        genetic_code_table: Genetic code table ID for translation (default: 1 = standard code)
    """

    def __init__(
        self,
        cds: dict[str, str],
        proteins: dict[str, str] | None = None,
        genetic_code_table: int | str = 1,
    ):
        """Initialize ReferenceSequenceSet with CDS and optional protein sequences.

        Args:
            cds: Dictionary mapping sequence IDs to CDS sequences
            proteins: Optional dictionary mapping sequence IDs to protein sequences
            genetic_code_table: Genetic code table ID (int) or name (str) for translation
        """
        self.cds = cds
        self.proteins = proteins
        self.genetic_code_table = genetic_code_table
        self._cds_strings_cache: list[str] | None = None
        self._protein_strings_cache: list[str] | None = None

    @classmethod
    def from_fasta(
        cls,
        cds_fasta_path: Path | str,
        protein_fasta_path: Path | str | None = None,
        genetic_code_table: int | str = 1,
    ) -> Self:
        """Construct ReferenceSequenceSet from FASTA file(s).

        Args:
            cds_fasta_path: Path to FASTA file containing CDS sequences
            protein_fasta_path: Optional path to FASTA file containing protein sequences
            genetic_code_table: Genetic code table ID (int) or name (str) for translation

        Returns:
            ReferenceSequenceSet instance

        Raises:
            FileNotFoundError: If FASTA file does not exist
            ValueError: If FASTA file is empty or contains duplicate IDs
        """
        # Load CDS sequences
        cds_path = Path(cds_fasta_path)
        if not cds_path.exists():
            raise FileNotFoundError(f"CDS FASTA file not found: {cds_path}")

        cds_dict: dict[str, str] = {}
        for record in SeqIO.parse(cds_path, "fasta"):
            if record.id in cds_dict:
                raise ValueError(f"Duplicate sequence ID in CDS FASTA: {record.id}")
            cds_dict[record.id] = str(record.seq)

        if not cds_dict:
            raise ValueError(f"CDS FASTA file is empty: {cds_path}")

        # Load protein sequences if provided
        proteins_dict: dict[str, str] | None = None
        if protein_fasta_path is not None:
            protein_path = Path(protein_fasta_path)
            if not protein_path.exists():
                raise FileNotFoundError(f"Protein FASTA file not found: {protein_path}")

            proteins_dict = {}
            for record in SeqIO.parse(protein_path, "fasta"):
                if record.id in proteins_dict:
                    raise ValueError(f"Duplicate sequence ID in protein FASTA: {record.id}")
                proteins_dict[record.id] = str(record.seq)

            if not proteins_dict:
                raise ValueError(f"Protein FASTA file is empty: {protein_path}")

        return cls(cds_dict, proteins_dict, genetic_code_table)

    def cds_strings(self, require_multiple_of_three: bool = True) -> list[str]:
        """Get normalized CDS strings.

        Normalization:
        - Convert to uppercase
        - Replace U with T (RNA to DNA)
        - Optionally validate length is multiple of 3

        Results are cached for efficiency.

        Args:
            require_multiple_of_three: If True, raise error if CDS length is not divisible by 3

        Returns:
            List of normalized CDS strings in consistent order

        Raises:
            ValueError: If require_multiple_of_three is True and any CDS length % 3 != 0
        """
        if self._cds_strings_cache is not None:
            # Return cached result, but still validate if required
            if require_multiple_of_three:
                self._validate_cds_lengths()
            return self._cds_strings_cache

        # Normalize and validate CDS sequences
        normalized: list[str] = []
        for seq_id, seq in self.cds.items():
            # Normalize: uppercase and U->T
            normalized_seq = seq.upper().replace("U", "T")

            # Validate length if required
            if require_multiple_of_three and len(normalized_seq) % 3 != 0:
                raise ValueError(
                    f"CDS sequence '{seq_id}' has length {len(normalized_seq)} "
                    f"which is not a multiple of 3 (index: {list(self.cds.keys()).index(seq_id)})"
                )

            normalized.append(normalized_seq)

        # Cache the result
        self._cds_strings_cache = normalized
        return normalized

    def _validate_cds_lengths(self) -> None:
        """Validate that all cached CDS sequences have length divisible by 3.

        Raises:
            ValueError: If any CDS length % 3 != 0
        """
        if self._cds_strings_cache is None:
            return

        for idx, (seq_id, cached_seq) in enumerate(
            zip(self.cds.keys(), self._cds_strings_cache, strict=True)
        ):
            if len(cached_seq) % 3 != 0:
                raise ValueError(
                    f"CDS sequence '{seq_id}' has length {len(cached_seq)} "
                    f"which is not a multiple of 3 (index: {idx})"
                )

    def protein_strings(
        self, strip_terminal_stop: bool = False, error_on_internal_stop: bool = True
    ) -> list[str]:
        """Get protein strings from provided proteins or by translating CDS.

        If proteins were provided during initialization, returns those.
        Otherwise, translates CDS sequences using Biopython with the specified genetic code.

        Results are cached for efficiency.

        Args:
            strip_terminal_stop: If True, remove terminal stop codon (*) from translations
            error_on_internal_stop: If True, raise error on internal stop codons

        Returns:
            List of protein strings in consistent order

        Raises:
            ValueError: If error_on_internal_stop is True and internal stop found
            ValueError: If CDS translation fails
        """
        if self._protein_strings_cache is not None:
            # Return cached result, but still validate if required
            if error_on_internal_stop:
                self._validate_no_internal_stops()
            return self._protein_strings_cache

        result: list[str] = []

        if self.proteins is not None:
            # Use provided protein sequences
            for seq_id in self.cds.keys():
                if seq_id not in self.proteins:
                    raise ValueError(
                        f"CDS sequence '{seq_id}' not found in provided protein sequences"
                    )
                protein_seq = self.proteins[seq_id]

                # Validate for internal stops if required
                if error_on_internal_stop:
                    stop_idx = protein_seq.find("*")
                    if stop_idx != -1 and stop_idx < len(protein_seq) - 1:
                        raise ValueError(
                            f"Protein sequence '{seq_id}' contains internal stop codon "
                            f"at position {stop_idx}"
                        )

                # Strip terminal stop if requested
                if strip_terminal_stop and protein_seq.endswith("*"):
                    protein_seq = protein_seq[:-1]

                result.append(protein_seq)
        else:
            # Translate CDS sequences
            cds_normalized = self.cds_strings(require_multiple_of_three=True)

            for idx, (seq_id, cds_seq) in enumerate(
                zip(self.cds.keys(), cds_normalized, strict=True)
            ):
                try:
                    # Translate using Biopython
                    bio_seq = Seq(cds_seq)
                    protein_seq = str(
                        bio_seq.translate(table=self.genetic_code_table)  # type: ignore[arg-type]
                    )

                    # Check for internal stops (stop codons before the last position)
                    if error_on_internal_stop:
                        stop_idx = protein_seq.find("*")
                        if stop_idx != -1 and stop_idx < len(protein_seq) - 1:
                            raise ValueError(
                                f"CDS sequence '{seq_id}' (index: {idx}) contains internal stop "
                                f"codon at protein position {stop_idx}"
                            )

                    # Strip terminal stop if requested
                    if strip_terminal_stop and protein_seq.endswith("*"):
                        protein_seq = protein_seq[:-1]

                    result.append(protein_seq)

                except Exception as e:
                    raise ValueError(
                        f"Failed to translate CDS sequence '{seq_id}' (index: {idx}): {e}"
                    ) from e

        # Cache the result
        self._protein_strings_cache = result
        return result

    def _validate_no_internal_stops(self) -> None:
        """Validate that cached protein sequences have no internal stop codons.

        Raises:
            ValueError: If any protein has internal stop codons
        """
        if self._protein_strings_cache is None:
            return

        for idx, (seq_id, protein_seq) in enumerate(
            zip(self.cds.keys(), self._protein_strings_cache, strict=True)
        ):
            stop_idx = protein_seq.find("*")
            if stop_idx != -1 and stop_idx < len(protein_seq) - 1:
                raise ValueError(
                    f"Protein sequence '{seq_id}' (index: {idx}) contains internal stop codon "
                    f"at position {stop_idx}"
                )

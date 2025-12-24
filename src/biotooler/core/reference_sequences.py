"""Reference sequence set for storing and processing CDS and protein sequences."""

from pathlib import Path
from typing import Self

import pandas as pd
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
        name: str | None = None,
    ):
        """Initialize ReferenceSequenceSet with CDS and optional protein sequences.

        Args:
            cds: Dictionary mapping sequence IDs to CDS sequences
            proteins: Optional dictionary mapping sequence IDs to protein sequences
            genetic_code_table: Genetic code table ID (int) or name (str) for translation
            name: Optional name identifier for this reference set
        """
        self.cds = cds
        self.proteins = proteins
        self.genetic_code_table = genetic_code_table
        self.name = name
        self._cds_strings_cache: list[str] | None = None
        self._cds_validated_multiple_of_three: bool = False
        self._protein_strings_cache: list[str] | None = None
        self._protein_validated_no_internal_stops: bool = False

    @classmethod
    def from_fasta(
        cls,
        cds_fasta_path: Path | str,
        protein_fasta_path: Path | str | None = None,
        genetic_code_table: int | str = 1,
        name: str | None = None,
    ) -> Self:
        """Construct ReferenceSequenceSet from FASTA file(s).

        Args:
            cds_fasta_path: Path to FASTA file containing CDS sequences
            protein_fasta_path: Optional path to FASTA file containing protein sequences
            genetic_code_table: Genetic code table ID (int) or name (str) for translation
            name: Optional name identifier for this reference set

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

        return cls(cds_dict, proteins_dict, genetic_code_table, name)

    @classmethod
    def from_csv(
        cls,
        csv_path: Path | str,
        id_column: str,
        cds_column: str,
        protein_column: str | None = None,
        genetic_code_table: int | str = 1,
        name: str | None = None,
    ) -> Self:
        """Construct ReferenceSequenceSet from CSV file.

        Args:
            csv_path: Path to CSV file
            id_column: Name of column containing sequence IDs
            cds_column: Name of column containing CDS sequences
            protein_column: Optional name of column containing protein sequences
            genetic_code_table: Genetic code table ID (int) or name (str) for translation
            name: Optional name identifier for this reference set

        Returns:
            ReferenceSequenceSet instance

        Raises:
            FileNotFoundError: If CSV file does not exist
            ValueError: If required columns are missing, CSV is empty, or contains duplicate IDs
        """
        # Load CSV file
        csv_file = Path(csv_path)
        if not csv_file.exists():
            raise FileNotFoundError(f"CSV file not found: {csv_file}")

        try:
            df = pd.read_csv(csv_file)
        except Exception as e:
            raise ValueError(f"Failed to read CSV file: {e}") from e

        if df.empty:
            raise ValueError(f"CSV file is empty: {csv_file}")

        # Validate required columns
        if id_column not in df.columns:
            raise ValueError(
                f"ID column '{id_column}' not found in CSV. "
                f"Available columns: {', '.join(df.columns)}"
            )
        if cds_column not in df.columns:
            raise ValueError(
                f"CDS column '{cds_column}' not found in CSV. "
                f"Available columns: {', '.join(df.columns)}"
            )

        # Check for protein column if specified
        if protein_column is not None and protein_column not in df.columns:
            raise ValueError(
                f"Protein column '{protein_column}' not found in CSV. "
                f"Available columns: {', '.join(df.columns)}"
            )

        # Extract CDS sequences
        cds_dict: dict[str, str] = {}
        # start=2 accounts for header row in 1-based row numbering
        for row_num, (_idx, row) in enumerate(df.iterrows(), start=2):
            seq_id = str(row[id_column])
            cds_value = row[cds_column]
            if pd.isna(cds_value):  # type: ignore[arg-type]
                raise ValueError(f"CDS sequence is missing for ID '{seq_id}' at row {row_num}")
            if seq_id in cds_dict:
                raise ValueError(f"Duplicate sequence ID in CSV: {seq_id}")
            cds_dict[seq_id] = str(cds_value)

        if not cds_dict:
            raise ValueError("No valid CDS sequences found in CSV")

        # Extract protein sequences if column is specified
        proteins_dict: dict[str, str] | None = None
        if protein_column is not None:
            proteins_dict = {}
            # start=2 accounts for header row in 1-based row numbering
            for row_num, (_idx, row) in enumerate(df.iterrows(), start=2):
                seq_id = str(row[id_column])
                protein_value = row[protein_column]
                if pd.isna(protein_value):  # type: ignore[arg-type]
                    raise ValueError(
                        f"Protein sequence is missing for ID '{seq_id}' at row {row_num}"
                    )
                proteins_dict[seq_id] = str(protein_value)

        return cls(cds_dict, proteins_dict, genetic_code_table, name)

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
            # Return cached result, but validate if required and not already validated
            if require_multiple_of_three and not self._cds_validated_multiple_of_three:
                self._validate_cds_lengths()
                self._cds_validated_multiple_of_three = True
            return self._cds_strings_cache

        # Normalize and validate CDS sequences
        normalized: list[str] = []
        for idx, (seq_id, seq) in enumerate(self.cds.items()):
            # Normalize: uppercase and U->T
            normalized_seq = seq.upper().replace("U", "T")

            # Validate length if required
            if require_multiple_of_three and len(normalized_seq) % 3 != 0:
                raise ValueError(
                    f"CDS sequence '{seq_id}' has length {len(normalized_seq)} "
                    f"which is not a multiple of 3 (index: {idx})"
                )

            normalized.append(normalized_seq)

        # Cache the result and mark validation state
        self._cds_strings_cache = normalized
        if require_multiple_of_three:
            self._cds_validated_multiple_of_three = True
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
            # Return cached result, but validate if required and not already validated
            if error_on_internal_stop and not self._protein_validated_no_internal_stops:
                self._validate_no_internal_stops()
                self._protein_validated_no_internal_stops = True
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

        # Cache the result and mark validation state
        self._protein_strings_cache = result
        if error_on_internal_stop:
            self._protein_validated_no_internal_stops = True
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

    def validate(self, kind: str = "both") -> None:
        """Validate sequences based on specified kind.

        Validates sequences without modifying the instance. Checks for empty sequences
        and kind-specific constraints (CDS length multiple of 3, no internal stops in proteins).

        Args:
            kind: Type of validation to perform:
                - "cds": Validate CDS sequences only (length multiple of 3, no empty)
                - "protein": Validate protein sequences only (no internal stops, no empty)
                - "both": Validate both CDS and protein sequences (default)

        Raises:
            ValueError: If kind is not one of "cds", "protein", or "both"
            ValueError: If any sequence is empty
            ValueError: If CDS validation fails (length not multiple of 3)
            ValueError: If protein validation fails (internal stop codons)
        """
        if kind not in ("cds", "protein", "both"):
            raise ValueError(f"Invalid kind '{kind}'. Must be one of: 'cds', 'protein', 'both'")

        # Validate CDS sequences
        if kind in ("cds", "both"):
            # Check for empty CDS sequences
            for seq_id, seq in self.cds.items():
                if not seq or not seq.strip():
                    raise ValueError(f"CDS sequence '{seq_id}' is empty")

            # Validate length is multiple of 3
            self.cds_strings(require_multiple_of_three=True)

        # Validate protein sequences
        if kind in ("protein", "both"):
            if self.proteins is not None:
                # Check for empty protein sequences
                for seq_id, seq in self.proteins.items():
                    if not seq or not seq.strip():
                        raise ValueError(f"Protein sequence '{seq_id}' is empty")

            # Validate no internal stops (this will use provided proteins or translate CDS)
            self.protein_strings(error_on_internal_stop=True)

    def with_proteins(self, proteins: dict[str, str]) -> Self:
        """Create a new ReferenceSequenceSet with updated protein sequences.

        Returns a new instance with the provided proteins merged with or replacing
        existing proteins. If a protein ID already exists with a different sequence,
        raises an error (no silent changes).

        Args:
            proteins: Dictionary mapping sequence IDs to protein sequences.
                Can include new proteins or replacements for existing ones.

        Returns:
            New ReferenceSequenceSet instance with updated proteins

        Raises:
            ValueError: If a protein ID exists with a different sequence
            ValueError: If any provided protein sequence is empty
        """
        # Check for empty sequences in provided proteins
        for seq_id, seq in proteins.items():
            if not seq or not seq.strip():
                raise ValueError(f"Provided protein sequence '{seq_id}' is empty")

        # Build new proteins dictionary
        new_proteins: dict[str, str] = {}
        if self.proteins is not None:
            new_proteins.update(self.proteins)

        # Check for conflicts before merging
        for seq_id, new_seq in proteins.items():
            if seq_id in new_proteins and new_proteins[seq_id] != new_seq:
                raise ValueError(
                    f"Protein sequence '{seq_id}' already exists with a different sequence. "
                    f"Existing: '{new_proteins[seq_id][:20]}...', "
                    f"Provided: '{new_seq[:20]}...'"
                )
            new_proteins[seq_id] = new_seq

        return self.__class__(self.cds, new_proteins, self.genetic_code_table)

    def extend_cds(self, cds: dict[str, str]) -> Self:
        """Create a new ReferenceSequenceSet with additional CDS sequences.

        Returns a new instance with the provided CDS sequences added. If a CDS ID
        already exists with a different sequence, raises an error (no silent changes).
        Cached data is invalidated in the new instance.

        Args:
            cds: Dictionary mapping sequence IDs to CDS sequences to add

        Returns:
            New ReferenceSequenceSet instance with extended CDS

        Raises:
            ValueError: If a CDS ID exists with a different sequence
            ValueError: If any provided CDS sequence is empty
        """
        # Check for empty sequences in provided CDS
        for seq_id, seq in cds.items():
            if not seq or not seq.strip():
                raise ValueError(f"Provided CDS sequence '{seq_id}' is empty")

        # Build new CDS dictionary
        new_cds: dict[str, str] = {}
        new_cds.update(self.cds)

        # Check for conflicts before extending
        for seq_id, new_seq in cds.items():
            if seq_id in new_cds:
                # Normalize both sequences for comparison
                existing_normalized = new_cds[seq_id].upper().replace("U", "T")
                new_normalized = new_seq.upper().replace("U", "T")
                if existing_normalized != new_normalized:
                    raise ValueError(
                        f"CDS sequence '{seq_id}' already exists with a different sequence. "
                        f"Existing: '{new_cds[seq_id][:20]}...', "
                        f"Provided: '{new_seq[:20]}...'"
                    )
                # If sequences match (after normalization), keep the existing one
            else:
                new_cds[seq_id] = new_seq

        return self.__class__(new_cds, self.proteins, self.genetic_code_table)

    def merge(self, other: Self) -> Self:
        """Merge this ReferenceSequenceSet with another.

        Returns a new instance containing sequences from both sets. Sequences from
        this instance take precedence in ordering. If the same ID exists in both
        sets with different sequences, raises an error (no silent changes).

        The genetic code table from this instance is preserved in the merged result.

        Args:
            other: Another ReferenceSequenceSet to merge with this one

        Returns:
            New ReferenceSequenceSet with merged sequences

        Raises:
            ValueError: If a sequence ID exists in both sets with different sequences
            ValueError: If any sequence in other is empty
        """
        # Check for empty sequences in other.cds
        for seq_id, seq in other.cds.items():
            if not seq or not seq.strip():
                raise ValueError(f"CDS sequence '{seq_id}' in other set is empty")

        # Check for empty sequences in other.proteins
        if other.proteins is not None:
            for seq_id, seq in other.proteins.items():
                if not seq or not seq.strip():
                    raise ValueError(f"Protein sequence '{seq_id}' in other set is empty")

        # Merge CDS sequences
        new_cds: dict[str, str] = {}
        new_cds.update(self.cds)

        for seq_id, other_seq in other.cds.items():
            if seq_id in new_cds:
                # Normalize both sequences for comparison
                existing_normalized = new_cds[seq_id].upper().replace("U", "T")
                other_normalized = other_seq.upper().replace("U", "T")
                if existing_normalized != other_normalized:
                    raise ValueError(
                        f"CDS sequence '{seq_id}' exists in both sets with different sequences. "
                        f"This: '{new_cds[seq_id][:20]}...', "
                        f"Other: '{other_seq[:20]}...'"
                    )
            else:
                new_cds[seq_id] = other_seq

        # Merge protein sequences if either set has proteins
        new_proteins: dict[str, str] | None = None
        if self.proteins is not None or other.proteins is not None:
            new_proteins = {}
            if self.proteins is not None:
                new_proteins.update(self.proteins)

            if other.proteins is not None:
                for seq_id, other_seq in other.proteins.items():
                    if seq_id in new_proteins:
                        if new_proteins[seq_id] != other_seq:
                            raise ValueError(
                                f"Protein sequence '{seq_id}' exists in both sets "
                                f"with different sequences. "
                                f"This: '{new_proteins[seq_id][:20]}...', "
                                f"Other: '{other_seq[:20]}...'"
                            )
                        # If sequences match, keep the existing one
                    else:
                        new_proteins[seq_id] = other_seq

        return self.__class__(new_cds, new_proteins, self.genetic_code_table)

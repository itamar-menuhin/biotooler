"""Error types for biotooler."""


class InvalidSequenceError(ValueError):
    """Raised when a sequence contains invalid characters for its molecule type.

    Attributes:
        message: Descriptive error message
        char: The invalid character that was found
        index: The position of the invalid character in the sequence
        molecule_type: The molecule type being validated (DNA, RNA, or protein)
    """

    def __init__(
        self,
        message: str,
        char: str | None = None,
        index: int | None = None,
        molecule_type: str | None = None,
    ):
        """Initialize InvalidSequenceError with detailed information.

        Args:
            message: Descriptive error message
            char: The invalid character that was found
            index: The position of the invalid character in the sequence
            molecule_type: The molecule type being validated
        """
        super().__init__(message)
        self.char = char
        self.index = index
        self.molecule_type = molecule_type

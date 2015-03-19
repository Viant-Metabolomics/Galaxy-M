""" Datatypes for Galaxy-M.
"""

from galaxy.datatypes.binary import (
    Binary,
    SQlite,
)


class SQliteSPS(SQlite):
    file_ext = "sps.sqlite"

Binary.register_sniffable_binary_format("sps.sqlite", "sps.sqlite", SQliteSPS)


class SQliteTM(SQlite):
    file_ext = "tm.sqlite"

Binary.register_sniffable_binary_format("tm.sqlite", "tm.sqlite", SQliteTM)


class SQliteEFS(SQlite):
    file_ext = "efs.sqlite"

Binary.register_sniffable_binary_format("efs.sqlite", "efs.sqlite", SQliteEFS)


class SQlitePPS(SQlite):
    file_ext = "pps.sqlite"

Binary.register_sniffable_binary_format("pps.sqlite", "pps.sqlite", SQlitePPS)


class SQliteOutput(SQlite):
    file_ext = "output.sqlite"

Binary.register_sniffable_binary_format("output.sqlite", "output.sqlite", SQliteOutput)

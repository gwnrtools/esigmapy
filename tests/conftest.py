import logging

import esigmapy.legacy as _legacy_module

# legacy.py uses logging.info without importing logging; patch it so tests work
_legacy_module.logging = logging

"""让 tests/ 下的用例能 import 上一层的模块（python/ 不是包，是扁平目录）。"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

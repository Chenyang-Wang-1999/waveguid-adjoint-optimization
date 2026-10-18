"""跑优化循环。

用法::

    cd data && python3 ../python/run_optimize.py

如果 ``data/DOF_list.csv`` 已存在，会从最后一轮继续（日志追加）。
"""

import optimizer
from paths import chdir_data


def main(**kwargs):
    chdir_data()
    return optimizer.run(**kwargs)


if __name__ == "__main__":
    main(iterations=100, alpha0=1e-3)

# ============================================================================
#  交叉波导伴随优化 —— 构建脚本
# ----------------------------------------------------------------------------
#  所有依赖路径都由 libMesh 自带的 libmesh-config 推导，仓库里不存绝对路径。
#  本机特有的设置写进 config.mk（不要提交，见 config.mk.example）。
#
#      make                 # 编译 C++ 求解器（默认 opt）
#      make METHOD=dbg      # 换 libMesh 的编译模式
#      make python          # 编译 Python 用的 Bezier 扩展模块（pybind11）
#      make run             # 在 data/ 里跑一次正问题
#      make docs            # 编译 Typst 文档
#      make clean
# ============================================================================

-include config.mk

# ---- 依赖位置 --------------------------------------------------------------
LIBMESH_DIR  ?= $(HOME)/local
LIBMESH_CONFIG := $(LIBMESH_DIR)/bin/libmesh-config

# libMesh 编译模式：opt / dbg / devel（取决于 libMesh 编译时的 configure 选项）
# 注意：libmesh-config 是从*环境变量*里读 METHOD 的，而 make 的 $(shell ...) 不会把
# makefile 变量（哪怕 export 过）传进去，所以下面每次调用都显式写成 METHOD=$(METHOD)。
METHOD ?= opt

# Python 解释器（编译 Bezier 扩展模块时用；带 pybind11 的那个）
PYTHON  ?= python3
PYBIND11_INC ?= $(shell $(PYTHON) -c "import pybind11; print(pybind11.get_include())" 2>/dev/null)
# 必须显式带上解释器自己的头文件目录：本机 /usr/include/python3.8 在编译器默认搜索路径里，
# 不加这一项会拿 3.8 的 Python.h 去编，虽然文件名后缀是 3.9，import 时仍会报版本不匹配。
PYTHON_INC ?= $(shell $(PYTHON)-config --includes 2>/dev/null)

# ---- 工具链 ----------------------------------------------------------------
ifeq ($(wildcard $(LIBMESH_CONFIG)),)
  $(error 找不到 $(LIBMESH_CONFIG)。请设置 LIBMESH_DIR，或在 config.mk 里写 LIBMESH_DIR=/你的/libmesh/安装前缀)
endif

CXX      := $(shell $(LIBMESH_CONFIG) --cxx)
CXXFLAGS := $(shell METHOD=$(METHOD) $(LIBMESH_CONFIG) --cxxflags) -std=gnu++17 -Wall
CPPFLAGS := $(shell METHOD=$(METHOD) $(LIBMESH_CONFIG) --include)
LIBS     := $(shell METHOD=$(METHOD) $(LIBMESH_CONFIG) --libs) -lmatio
INCLUDES := -I$(CURDIR)/waveguide_optimization

# ---- 目标 ------------------------------------------------------------------
SRCDIR   := waveguide_optimization
SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS  := $(SOURCES:.cpp=.o)
DEPS     := $(OBJECTS:.o=.d)
TARGET   := $(SRCDIR)/waveguide_optimization

# Python 扩展模块（Bezier.cpp 里的 pybind11 绑定，供 python/original_mesh.py 使用）
PYEXT    := $(shell $(PYTHON)-config --extension-suffix 2>/dev/null || echo .so)
PYTARGET := python/Bezier$(PYEXT)

.PHONY: all run python docs clean distclean help

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $^ -o $@ $(LIBS)

$(SRCDIR)/%.o: $(SRCDIR)/%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

# Bezier 扩展模块：只编 Bezier.cpp，不链接 libMesh
python: $(PYTARGET)

$(PYTARGET): $(SRCDIR)/Bezier.cpp $(SRCDIR)/Bezier.h $(SRCDIR)/TriDiagonal.h
	@test -n "$(PYBIND11_INC)" || (echo "找不到 pybind11，请先 pip install pybind11" && false)
	$(CXX) -O2 -std=gnu++17 -shared -fPIC -DEXPORT_PYLIB \
	    $(PYTHON_INC) -I$(PYBIND11_INC) $(INCLUDES) \
	    -I$(LIBMESH_DIR)/eigen \
	    $(SRCDIR)/Bezier.cpp -o $@

# 跑一次正问题（输入输出都在 data/，输出文件已被 .gitignore 忽略）
run: $(TARGET)
	cd data && mpiexec -n 1 ../$(TARGET) $(ARGS)

docs:
	./docs/build.sh

clean:
	rm -f $(OBJECTS) $(DEPS) $(TARGET) $(PYTARGET)
	rm -rf __pycache__ python/__pycache__

distclean: clean
	rm -rf docs/build
	rm -f data/out_*.msh data/*_result.mat

help:
	@echo "make                 编译求解器 (METHOD=opt)"
	@echo "make METHOD=dbg      用 libMesh 的 dbg 模式编译"
	@echo "make python          编译 python/Bezier*.so"
	@echo "make run ARGS=...    在 data/ 里跑一次"
	@echo "make docs            编译 Typst 文档"
	@echo "make clean           清理编译产物"
	@echo ""
	@echo "LIBMESH_DIR = $(LIBMESH_DIR)"
	@echo "CXX         = $(CXX)"
	@echo "METHOD      = $(METHOD)"

-include $(DEPS)

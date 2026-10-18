// HDF5 读写小工具：给 C++ 后端与 Python 前端之间传数据用。
//
// 设计约定（两边必须一致）：
//   * 实数矩阵    -> 二维数据集 (rows, cols)，FLOAT64
//   * 整数矩阵    -> 二维数据集 (rows, cols)，HDF5 读时会自动转换类型，
//                    所以 Python 写 int32/int64 都能读
//   * 复数向量    -> 二维数据集 (n, 2)，FLOAT64，实部在 [:,0]、虚部在 [:,1]，
//                    并带属性 complex = 1。这与 PETSc 的 HDF5 viewer 写复数
//                    Vec 的约定一致，所以将来若改用 VecView 写场，Python 侧
//                    的读取代码不用变。
//   * 标量/字符串 -> 组的属性（attributes），放在 /meta 组里
//
// 依赖：hdf5.h + libhdf5。两者本来就在 libmesh-config 的 --include / --libs
// 里，所以不需要额外改编译命令。
//
// 错误处理：HDF5 默认静默返回负值，这里统一转成异常。

#ifndef H5IO_H_INCLUDED
#define H5IO_H_INCLUDED

#include <hdf5.h>

#include <Eigen/Core>
#include <complex>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace h5io
{

// result.h5 的数据布局版本。任何会让旧文件读不出来的改动都必须把它 +1。
//   1: 每 rank 一份（S_params_<rank> / Ez_<rank>）
//   2: 各 rank 在 C++ 侧 MPI 归约后只写一份
static const int FORMAT_VERSION = 2;

inline void fail(const std::string & what)
{
  throw std::runtime_error("h5io: " + what);
}

inline hid_t ok(hid_t id, const std::string & what)
{
  if (id < 0)
    fail(what);
  return id;
}

inline void ok(herr_t status, const std::string & what)
{
  if (status < 0)
    fail(what);
}

// 关掉 HDF5 自带的报错打印，改由异常报告
inline void quiet_errors()
{
  H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
}

inline bool file_exists(const std::string & name)
{
  if (name.empty())
    return false;
  std::ifstream f(name);
  return f.good();
}

// RAII 文件句柄
class File
{
public:
  explicit File(const std::string & name, bool for_write)
  {
    quiet_errors();
    _fid = ok(for_write
                ? H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)
                : H5Fopen(name.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT),
              "open " + name);
  }

  // 打开已存在的文件用于追加写（MPI 下按 rank 轮流写时用这个）
  struct append_t {};
  File(const std::string & name, append_t)
  {
    quiet_errors();
    _fid = ok(H5Fopen(name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT), "open(rw) " + name);
  }

  ~File()
  {
    if (_fid >= 0)
      H5Fclose(_fid);
  }

  File(const File &) = delete;
  File & operator=(const File &) = delete;

  hid_t id() const { return _fid; }

private:
  hid_t _fid = -1;
};

// 取一个组，不存在就创建
inline hid_t group(hid_t loc, const std::string & path)
{
  htri_t has = H5Lexists(loc, path.c_str(), H5P_DEFAULT);
  if (has > 0)
    return ok(H5Gopen2(loc, path.c_str(), H5P_DEFAULT), "open group " + path);
  return ok(H5Gcreate2(loc, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
            "create group " + path);
}

class Group
{
public:
  Group(hid_t loc, const std::string & path) : _gid(group(loc, path)) {}
  ~Group()
  {
    if (_gid >= 0)
      H5Gclose(_gid);
  }
  Group(const Group &) = delete;
  Group & operator=(const Group &) = delete;
  hid_t id() const { return _gid; }

private:
  hid_t _gid = -1;
};

inline void unlink_if_exists(hid_t loc, const std::string & name)
{
  H5E_BEGIN_TRY
  {
    H5Ldelete(loc, name.c_str(), H5P_DEFAULT);
  }
  H5E_END_TRY
}

// ---------------------------------------------------------------- 写入

inline void write_matrix(hid_t loc, const std::string & name, const Eigen::MatrixXd & m)
{
  unlink_if_exists(loc, name);
  hsize_t dims[2] = {(hsize_t)m.rows(), (hsize_t)m.cols()};
  hid_t sp = ok(H5Screate_simple(2, dims, nullptr), "dataspace " + name);
  hid_t ds = ok(H5Dcreate2(loc, name.c_str(), H5T_NATIVE_DOUBLE, sp,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                "create " + name);
  // Eigen 默认列主序，这里按行主序展开写，Python 侧读到的就是 (rows, cols)
  std::vector<double> buf((std::size_t)m.rows() * m.cols());
  for (Eigen::Index i = 0; i < m.rows(); ++i)
    for (Eigen::Index j = 0; j < m.cols(); ++j)
      buf[(std::size_t)i * m.cols() + j] = m(i, j);
  ok(H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()),
     "write " + name);
  H5Dclose(ds);
  H5Sclose(sp);
}

inline void write_matrix(hid_t loc, const std::string & name, const Eigen::MatrixXi & m)
{
  unlink_if_exists(loc, name);
  hsize_t dims[2] = {(hsize_t)m.rows(), (hsize_t)m.cols()};
  hid_t sp = ok(H5Screate_simple(2, dims, nullptr), "dataspace " + name);
  hid_t ds = ok(H5Dcreate2(loc, name.c_str(), H5T_NATIVE_INT, sp,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                "create " + name);
  std::vector<int> buf((std::size_t)m.rows() * m.cols());
  for (Eigen::Index i = 0; i < m.rows(); ++i)
    for (Eigen::Index j = 0; j < m.cols(); ++j)
      buf[(std::size_t)i * m.cols() + j] = m(i, j);
  ok(H5Dwrite(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()),
     "write " + name);
  H5Dclose(ds);
  H5Sclose(sp);
}

inline void set_attr(hid_t loc, const std::string & name, int value)
{
  unlink_if_exists(loc, name);
  hid_t sp = ok(H5Screate(H5S_SCALAR), "attr dataspace " + name);
  hid_t at = ok(H5Acreate2(loc, name.c_str(), H5T_NATIVE_INT, sp,
                           H5P_DEFAULT, H5P_DEFAULT),
                "attr " + name);
  ok(H5Awrite(at, H5T_NATIVE_INT, &value), "write attr " + name);
  H5Aclose(at);
  H5Sclose(sp);
}

inline void set_attr(hid_t loc, const std::string & name, double value)
{
  unlink_if_exists(loc, name);
  hid_t sp = ok(H5Screate(H5S_SCALAR), "attr dataspace " + name);
  hid_t at = ok(H5Acreate2(loc, name.c_str(), H5T_NATIVE_DOUBLE, sp,
                           H5P_DEFAULT, H5P_DEFAULT),
                "attr " + name);
  ok(H5Awrite(at, H5T_NATIVE_DOUBLE, &value), "write attr " + name);
  H5Aclose(at);
  H5Sclose(sp);
}

inline void set_attr(hid_t loc, const std::string & name, const std::string & value)
{
  unlink_if_exists(loc, name);
  hid_t st = ok(H5Tcopy(H5T_C_S1), "string type");
  ok(H5Tset_size(st, H5T_VARIABLE), "string size");
  hid_t sp = ok(H5Screate(H5S_SCALAR), "attr dataspace " + name);
  hid_t at = ok(H5Acreate2(loc, name.c_str(), st, sp, H5P_DEFAULT, H5P_DEFAULT),
                "attr " + name);
  const char * c = value.c_str();
  ok(H5Awrite(at, st, &c), "write attr " + name);
  H5Aclose(at);
  H5Sclose(sp);
  H5Tclose(st);
}

// 复数向量 -> (n,2) 实数 + 属性 complex=1
inline void write_complex_vector(hid_t loc, const std::string & name,
                                 const Eigen::VectorXcd & v)
{
  unlink_if_exists(loc, name);
  hsize_t dims[2] = {(hsize_t)v.size(), 2};
  hid_t sp = ok(H5Screate_simple(2, dims, nullptr), "dataspace " + name);
  hid_t ds = ok(H5Dcreate2(loc, name.c_str(), H5T_NATIVE_DOUBLE, sp,
                           H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                "create " + name);
  std::vector<double> buf(2 * (std::size_t)v.size());
  for (Eigen::Index i = 0; i < v.size(); ++i)
    {
      buf[2 * (std::size_t)i] = v(i).real();
      buf[2 * (std::size_t)i + 1] = v(i).imag();
    }
  ok(H5Dwrite(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()),
     "write " + name);
  set_attr(ds, "complex", 1);
  H5Dclose(ds);
  H5Sclose(sp);
}

// ---------------------------------------------------------------- 读取

namespace detail
{

// 读出数据集的形状，并返回一个已打开的 dataspace
inline hid_t shape_of(hid_t loc, const std::string & name, hsize_t dims[2])
{
  hid_t ds = ok(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), "open dataset " + name);
  hid_t sp = ok(H5Dget_space(ds), "get space " + name);
  int nd = H5Sget_simple_extent_ndims(sp);
  if (nd == 2)
    ok(H5Sget_simple_extent_dims(sp, dims, nullptr), "dims " + name);
  else if (nd == 1)
    {
      hsize_t d = 0;
      ok(H5Sget_simple_extent_dims(sp, &d, nullptr), "dims " + name);
      dims[0] = 1;
      dims[1] = d;
    }
  else
    fail("unsupported rank for " + name);
  H5Dclose(ds);
  return sp;
}

} // namespace detail

inline Eigen::MatrixXd read_matrix(hid_t loc, const std::string & name)
{
  hsize_t dims[2] = {0, 0};
  hid_t sp = detail::shape_of(loc, name, dims);
  hid_t ds = ok(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), "reopen " + name);
  std::vector<double> buf((std::size_t)dims[0] * dims[1]);
  ok(H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()),
     "read " + name);
  Eigen::MatrixXd m((Eigen::Index)dims[0], (Eigen::Index)dims[1]);
  for (hsize_t i = 0; i < dims[0]; ++i)
    for (hsize_t j = 0; j < dims[1]; ++j)
      m((Eigen::Index)i, (Eigen::Index)j) = buf[(std::size_t)i * dims[1] + j];
  H5Dclose(ds);
  H5Sclose(sp);
  return m;
}

inline Eigen::MatrixXi read_matrix_int(hid_t loc, const std::string & name)
{
  hsize_t dims[2] = {0, 0};
  hid_t sp = detail::shape_of(loc, name, dims);
  hid_t ds = ok(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), "reopen " + name);
  std::vector<int> buf((std::size_t)dims[0] * dims[1]);
  // HDF5 会自动把存储类型（int32/int64）转换成这里请求的 NATIVE_INT
  ok(H5Dread(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()),
     "read " + name);
  Eigen::MatrixXi m((Eigen::Index)dims[0], (Eigen::Index)dims[1]);
  for (hsize_t i = 0; i < dims[0]; ++i)
    for (hsize_t j = 0; j < dims[1]; ++j)
      m((Eigen::Index)i, (Eigen::Index)j) = buf[(std::size_t)i * dims[1] + j];
  H5Dclose(ds);
  H5Sclose(sp);
  return m;
}

// (n,2) 实数 + complex=1 -> 复数向量
inline Eigen::VectorXcd read_complex_vector(hid_t loc, const std::string & name)
{
  hsize_t dims[2] = {0, 0};
  hid_t sp = detail::shape_of(loc, name, dims);
  hid_t ds = ok(H5Dopen2(loc, name.c_str(), H5P_DEFAULT), "reopen " + name);
  std::vector<double> buf((std::size_t)dims[0] * dims[1]);
  ok(H5Dread(ds, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data()),
     "read " + name);
  Eigen::Index n = (Eigen::Index)(dims[1] == 2 ? dims[0] : dims[0] * dims[1]);
  Eigen::VectorXcd v(n);
  for (Eigen::Index i = 0; i < n; ++i)
    v(i) = dims[1] == 2
             ? std::complex<double>(buf[2 * (std::size_t)i], buf[2 * (std::size_t)i + 1])
             : std::complex<double>(buf[(std::size_t)i], 0.0);
  H5Dclose(ds);
  H5Sclose(sp);
  return v;
}

inline bool has(hid_t loc, const std::string & name)
{
  return H5Lexists(loc, name.c_str(), H5P_DEFAULT) > 0;
}

} // namespace h5io

#endif // H5IO_H_INCLUDED

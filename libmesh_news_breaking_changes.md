# libMesh NEWS digest: client-code-breaking changes

## Version range

`1.6.2 -> 1.7.0` through `1.8.1 -> 1.8.2` (i.e. libMesh 1.6.2 → 1.9.0-pre/devel); NEWS lines 1–612 of `/home/wangchenyang/install/libmesh/NEWS`, read from the "1.8.1 -> 1.8.2" heading at the top down through the `1.6.2 -> 1.7.0` section.

Method note: every bullet below is anchored to a verbatim NEWS line (quoted in *Quotes*, with line numbers). Where I also confirmed the claim against the checked-out sources (`devel`, and `origin/branch_v1.7.x` / `origin/branch_v1.8.x`) I say so explicitly as **[verified in header]**. Anything I could not confirm is listed as unconfirmed rather than asserted.

---

## Breaking changes that matter to client code

- **C++17 is now mandatory.** Under `1.7.3 -> 1.8.0`, Misc: "C++17 support is now required". *Action:* compile client code with `-std=c++17` (or later). No explicit minimum compiler version is stated anywhere in this range — only the C++17 requirement, plus the separate "Require support for enum forward declarations" (also `1.7.3 -> 1.8.0`).

- **`MeshTools::find_boundary_nodes()` overload removed.** Under `1.7.3 -> 1.8.0`: "Remove long-deprecated overload `MeshTools::find_boundary_nodes()`". *Action:* switch to the surviving overload; code calling the removed one no longer compiles.

- **`libmesh_make_unique` deprecated, `UniquePtr` deprecated, and the deprecated `UniquePtr` workaround removed.** Under `1.7.3 -> 1.8.0`: "Use `std::make_unique` everywhere, deprecate `libmesh_make_unique`"; "Remove deprected `UniquePtr` workaround for `std::unique_ptr`" [sic, "deprected" in the original]; and under `1.6.2 -> 1.7.0`: "Officially deprecate `UniquePtr`". *Action:* use `std::unique_ptr` / `std::make_unique`; expect `UniquePtr` to disappear.

- **`Elem::centroid()` deprecated in favour of new names.** Under `1.6.2 -> 1.7.0`: "Add `Elem::true_centroid()`, `Elem::vertex_average()`, deprecate `Elem::centroid()`". *Action:* replace `centroid()` calls; the deprecation implies removal later (behavioural note: `true_centroid()` and `vertex_average()` are *not* the same quantity).

- **`Side` and `SideEdge` classes deprecated; proxy sides are no longer built by default.** Under `1.6.2 -> 1.7.0`: "Deprecate the `SideEdge` class", "Deprecate the `Side` class", "Deprecate attempts to build proxy `Side` objects", "Stop building proxy sides by default". *Action:* client code that obtained proxy side/edge objects from `Elem` must move to the `build_side_ptr`/`build_edge_ptr` path (`"Add const version of build_edge_ptr"`, `"Add \"fill\" overload for build_edge_ptr"` are also in that section). This is the change most likely to silently change what a pointer points to.

- **Public `MeshBase::boundary_info` member is deprecated.** Under `1.6.2 -> 1.7.0`: "Deprecate public `MeshBase::boundary_info` member", plus "Add global data structures for mesh subdomain and boundary ids". **[verified in header]** `mesh_base.h` now declares `std::unique_ptr<BoundaryInfo> boundary_info;` with the comment *"Direct access to this class is now officially deprecated and will be removed in future libMesh versions. Use the `get_boundary_info()` accessor instead."* *Action:* migrate `mesh.boundary_info...` to `mesh.get_boundary_info()...`. Note the member being a `unique_ptr` means it is already `mesh.boundary_info->...`, not `mesh.boundary_info....` — but that pointer form predates this range (**[verified]** it is already `std::unique_ptr<BoundaryInfo>` at `branch_v1.7.x`), so the `.`→`->` adjustment is *not* a 1.6.2→1.9.0 break; only the deprecation is.

- **`BoundaryInfo` id handling gained new APIs and a new flag.** Under `1.7.3 -> 1.8.0`: "Add `BoundaryInfo::renumber_id()`" and "Add option for removing global boundary ID in `remove_id()`". **[verified in header]** `void remove_id (boundary_id_type id, bool global = false);` and `void renumber_id (boundary_id_type old_id, boundary_id_type new_id);`. *Action:* existing `remove_id(id)` calls still compile via the new default argument, but the default is `global = false`. I **could not confirm** whether `false` reproduces pre-1.8.0 behaviour or whether callers who previously relied on the global id-to-name entry being deleted must now pass `global = true` — check `remove_id()` call sites if you use boundary id *names*.

- **`exodusII.h` is no longer pulled in by libMesh headers.** Under `1.6.2 -> 1.7.0`, ExodusIO: "Get `exodusII.h` include out of our headers". *Action:* client code that used Exodus types only because a libMesh header transitively included `exodusII.h` must now `#include "exodusII.h"` (or avoid the type) — a real compile break with no deprecation cycle.

- **`ExodusII_IO_Helper::nodal_var_values` changed container type.** Under `1.6.2 -> 1.7.0`: "Make `ExodusII_IO_Helper::nodal_var_values` a map". *Action:* client code indexing it as a vector no longer compiles; rework to map access. (Also in that section: "`inquire()` function signature change, move to anonymous namespace" — the NEWS entry does not name the class or the new signature, so I **cannot confirm** the client-visible impact; it reads as an internal helper change.)

- **`NumericVector::restore_subvector()` takes a `std::unique_ptr` by value.** Under `1.7.3 -> 1.8.0`: "Pass `unique_ptr` by value to `NumericVector::restore_subvector()`". **[verified in header]** `virtual void restore_subvector(std::unique_ptr<NumericVector<T>>, const std::vector<numeric_index_type> &)`. *Action:* any `NumericVector<T>` subclass that overrides `restore_subvector` must change its signature and add `override` (without `override` the old raw-pointer signature silently stops overriding). **[verified]** this API does not exist at `branch_v1.7.x`, so it is new-in-range rather than a break for ≤1.7 code.

- **`NumericVector` gained a new *pure virtual* `max_allowed_id()`.** Under `1.6.2 -> 1.7.0`: "Add virtual `NumericVector::max_allowed_id()` API." **[verified in header]** `virtual std::size_t max_allowed_id() const = 0;`. *Action:* every custom `NumericVector<T>` subclass compiled against ≤1.6.2 must implement it or it cannot be instantiated. (Confirmed absent → present: it appears first in the 1.7 line, present at `branch_v1.7.x`.)

- **`linear_solver` member moved up the hierarchy.** Under `1.6.2 -> 1.7.0`: "Move `linear_solver` member from `LinearImplicitSystem` to `ImplicitSystem`". **[verified in header]** `mutable std::unique_ptr<LinearSolver<Number>> linear_solver;` is declared in `implicit_system.h`'s `ImplicitSystem`. *Action:* accessors keep working through inheritance, but (a) code in a `LinearImplicitSystem` subclass that declared its own `linear_solver` now conflicts/shadow-warns, and (b) this is an ABI break for pre-built client objects.

- **`ImplicitSystem::get_linear_solver()` reportedly now throws.** Under `1.6.2 -> 1.7.0`: "`ImplicitSystem::get_linear_solver()` now throws an error". *Action if true:* client code that tested the returned pointer for null and fell back must switch to exception handling. **Unconfirmed:** I read both implementations in the 1.9.0-pre sources and found no `libmesh_error`/`libmesh_not_implemented` in `ImplicitSystem::get_linear_solver()` (`src/systems/implicit_system.C`) — it resets, rebuilds via `LinearSolver<Number>::build()`, calls `init_systems(*this)`, and returns the pointer — and `LinearImplicitSystem::get_linear_solver()` still just returns `linear_solver.get()`. Treat the NEWS line as the statement of intent and verify against your libMesh build.

- **`System::is_initialized()` became `const`.** Under `1.7.3 -> 1.8.0`: "Make `System::is_initialized()` const". *Action:* any client subclass that overrides it must add `const` (a non-const override stops overriding — silently, unless `override` is used).

- **`SparsityPattern::Build` direct access deprecated, and `Build::mesh` removed.** Under `1.6.2 -> 1.7.0`: "Deprecate `SparsityPattern::Build` direct access", "Remove `SparsityPattern::Build::mesh`", "default+delete `SparsityPattern::Build` functions". *Action:* client code that reached into `SparsityPattern::Build` members must use the new accessors (`"Add SparsityPattern::build accessors, reorganize"`, "`DofMap::get_sparsity_pattern()` accessor", "Make `DofMap::build_sparsity()` public" — all same section). "Use SparsityPattern for `n_nz` in `SparseMatrix`" is in the same section, but **[verified]** `sparse_matrix.h` exposes no `n_nz`/`n_oz` accessor in either `branch_v1.7.x` or `devel`, so I found **no client-visible signature break** there.

- **`CouplingMatrix` index type changed.** Under `1.7.3 -> 1.8.0`: "Change `unsigned int` to `size_t` in `CouplingMatrix` to prevent overflow for large number of variables (#3541)". *Action:* client code storing/looping coupling indices in `unsigned int` still compiles but now narrows (warnings, and silent truncation in huge-variable-count runs); code taking `&CouplingMatrix::operator()` or relying on the exact signature breaks.

- **Ownership/type changes in several public members and ctor arguments** (all `1.6.2 -> 1.7.0`): "`UnsteadySolver: old_local_nonlinear_solution` is now a `std::shared_ptr`", "`PeriodicBoundaries` now stores `unique_ptrs` instead of dumb ptrs", "`MeshFunction: use unique_ptr` to manage `PointLocator`", "`PointLocatorTree: use std::shared_ptr` to manage tree", "`ParsedFEMFunction`: change \"parsers\" to a vector of unique_ptrs", "`System: store unique_ptrs in _vectors` map". *Action:* client code touching these members (including derived-class code iterating `_vectors`) must update types; any code that copied these containers can no longer do so.

- **`PetscMatrix<T>` is no longer the root of its own hierarchy.** Under `1.7.3 -> 1.8.0`: "Add `PetscMatrixBase<T>`, base class for `PetscMatrix<T>` and `PetscMatrixShellMatrix<T>`" (also "Add `StaticCondensation` class, derived from `PetscMatrixShellMatrix<T>`"). *Action:* `dynamic_cast`/`typeid` logic and any `template` code specialized strictly on `PetscMatrix<T>` may need to move to `PetscMatrixBase<T>`; partial specializations on the old base can silently stop matching. Related, area-11 note: a shell-matrix type now exists (`PetscMatrixShellMatrix<T>`), and "Add shell matrix support for `CondensedEigenSystem`".

- **Reduced-Basis API breaks (outside your 12 areas, listing for completeness).** Under `1.7.3 -> 1.8.0`: "`RBParameters`: store `std::vector<Real>` instead of single `Real`", "Deprecate `RBParameters::get_parameters_map()`", "Replace `NumericVector` with `std::vector` in `RBConstructionBase`", "const correctness updates for `RBThetaExpansion`, `RBConstruction` APIs". *Action:* RB client code must be recompiled and re-checked; a parameter that used to be a scalar is now a vector.

---

## Behavioural (silent) changes

These still compile; results or runtime behaviour can change.

- **`DofMap::distribute_dofs()` now returns the *global* DOF count.** `1.6.2 -> 1.7.0`: "`DofMap::distribute_dofs()` now returns the total number of DOFs across all procs". **[verified in header]** the return type is `std::size_t` both before and after, so **this compiles silently** and only the value changes — in parallel, code doing `n = dof_map.distribute_dofs(mesh)` now gets the global total where it previously got something local. Highest-risk silent change in the range.

- **`FEMContext` variable iteration order changed.** `1.6.2 -> 1.7.0`, FEMContext group: "Sort `FEMContext::_active_vars`", "Use `_active_vars` in FEMContext calculations", "Optional active_vars in all `FEMContext::FEMContext`", "Add active_vars member, accessor to FEMContext". Variables are now sorted and actually used to restrict calculations, so code that assumed the caller's ordering (hand-written index math, coupling matrices sized by position) can silently compute different quantities. **[verified in header]** ctors take `const std::vector<unsigned int> * active_vars = nullptr`, accessor `active_vars()`; `get_element_fe`, `get_side_fe`, `elem_fe_reinit(const std::vector<Point> * const pts = nullptr)`, `side_fe_reinit()` and `get_element_qrule()` (returning `const QBase &`) are all still present with these shapes.

- **`ASPECT_RATIO` was re-purposed; the old generic metric was renamed.** `1.7.3 -> 1.8.0`: "Rename existing, generic `ASPECT_RATIO` metric to `EDGE_LENGTH_RATIO`" plus "Add improved `ASPECT_RATIO` quality metric for Quad elements with unit tests". **[verified in header]** `ASPECT_RATIO` exists at `branch_v1.7.x` while `EDGE_LENGTH_RATIO` does not, and both exist in `devel`. So `elem->quality(ASPECT_RATIO)` still compiles and can return different numbers than before. *Action:* if you wanted the old generic edge-length ratio, switch to `EDGE_LENGTH_RATIO`.

- **`subdomain_ids()` semantics changed.** `1.6.2 -> 1.7.0`: "Include unpartitioned elements in `MeshBase::subdomain_ids`".

- **Dirichlet BC consistency checking is on again.** `1.6.2 -> 1.7.0`: "Re-enable `check_dirichlet_bcid_consistency`" — previously-silent inconsistent BC ids can now trip a check/error.

- **Boundary info now also lives on child elements.** `1.7.3 -> 1.8.0`: "Allow child elements to have `BoundaryInfo` associated with them". Boundary queries/iteration can now surface entries on child (refined) elements that were previously invisible.

- **Side/edge elements now carry `subdomain_id` and `p_level`.** `1.6.2 -> 1.7.0`: "Set `subdomain_id+p_level` in all side/edge elems" — code reading `subdomain_id()`/`p_level()` off side elems sees new values.

- **Boundary id-to-name maps are now consulted when regenerating id sets.** `1.6.2 -> 1.7.0`, BoundaryInfo: "Act on id to name maps in `BoundaryInfo::regenerate_id_sets`", and "Selectively remove boundary information after mesh stitching", "Add set names from other mesh when stitching". Relevant to `get_id_by_name()` results after `regenerate_id_sets()` / stitching.

- **`NumericVector` indexing add/set re-opens the vector.** `1.6.2 -> 1.7.0`: "`NumericVector`: Make indexing add/set operations toggle `_is_closed`" — code that cached a closed state or relied on no implicit re-assembly/close may now pay for (or observe) a re-close. Same section: "DistributedVector `GHOSTED` in serial can be `SERIAL`" — the type you get in a `--disable-mpi` build can differ.

- **`EigenSystem::solve` now errors when constraints are present.** `1.7.3 -> 1.8.0`: "Throw an error in `EigenSystem::solve` if we have constraints" — runs that used to proceed (possibly wrongly) now abort.

- **`(Equation)System(s)::init()` is re-implemented on top of `reinit_mesh()`.** `1.7.3 -> 1.8.0`: "Implement `(Equation)System(s)::init()` in terms of `reinit_mesh()`". **[verified in header]** `System::reinit_mesh()` is a public `virtual`. Silent risk: mesh-reinit side effects now occur on the `init()` path, and a client subclass that already had its own `reinit_mesh()` now overrides a library virtual.

- **Mesh/DOF bookkeeping fixes that change numbers.** `1.7.3 -> 1.8.0`: "Update id counts in `MeshCommunication::broadcast()`", "Fix `ReplicatedMesh(DistributedMesh&)` construction"; `1.8.1 -> 1.8.2`: "Don't look for 2D sides of a 2D element" (side iteration on 2D elements changes — directly relevant to `BoundaryInfo`/side loops), "Fix quadratic element issue in `TriangulatorInterface`", "Don't let GMSH physicals overwrite each other (#4128)"; `1.7.3 -> 1.8.0`: "Fix `Elem::is_child_on_edge()`", "Error if mismatched p-levels for Lagrange between neighbors" (new hard error), "Average side nodal solutions for H_DIV elements", "Remove caching from `FE<CLOUGH>`", "Fix `NumericVector::compare()`".

- **Thread-safety changes can change timing/behaviour but not results.** `1.8.1 -> 1.8.2`: "Lock at `PetscMatrix` level when doing add/set (#4113)", "Add mutex to `setEpsilon` (#4111)", "Fix double-creation of `PetscViewers`".

---

## Things I checked that are NOT breaking

- **`Parameters::get<T>()` return type — unchanged, and the `std::set` line does not affect it.** This is the one you flagged, so, precisely:

  - `get<T>()` returns **`const T &`** in *both* the old and new trees: **[verified]** `template <typename T> const T & get (std::string_view) const;` in `devel`, and the identical `const T & get (const std::string &) const;` at `origin/branch_v1.7.x`. It was never by value and never a non-const reference, and it is documented as *"A constant reference to the specified parameter value"*. `set<T>()` returns `T &` (*"A writable reference"*), `insert<T>()` returns `void` (*"does not return a writable reference"*).
  - **The NEWS in this entire range contains no line about `get()` returning by value, about a `const T &`, about dangling references, or about changing reference-ness.** "Allow std::set Parameter instances." (`1.8.1 -> 1.8.2`) is about which *value types* a parameter may hold/print: **[verified]** the only related header delta is a `print_helper(std::ostream &, const std::set<P, C, A> * param)` overload (present at `branch_v1.8.x` and `devel`), i.e. a `std::set` parameter is now printable instead of failing to compile in `Parameters::print()`. It does not touch `get`/`set`/`have_parameter` signatures.
  - **What that implies for `double & p = params.get<double>("x");`** — this is a **compile error, then and now**: `get<double>` yields `const double &`, and a non-const `double &` cannot bind to it. It is not, and never was, a dangling-reference bug; it simply does not compile. Correct forms: `const double & p = params.get<double>("x");` or `double p = params.get<double>("x");`.
  - **The real dangling-reference hazard is elsewhere** (and is *not* mentioned in the NEWS): **[verified in header]** `Parameters::set<T>()` stores into a `std::map` of `std::unique_ptr<ParameterBase>` and does `_values[name] = std::make_unique<Parameter<T>>()` whenever the existing entry is not already a `Parameter<T>`. So calling `set<T>()` on an existing name with a *different* type destroys the old `Parameter<T>` object — any `const T &` you obtained earlier from `get<T>()`/`set<T>()` for that name then dangles. `remove(name)` (`_values.erase`) and `operator+=` (documented: *"Any parameters of the same name already in \p this are replaced"*) have the same effect. Hold references only across reads, or re-`get` after any mutation.
- **"Avoid out of bounds access w/ operator[] and assert correctness."** (`1.8.1 -> 1.8.2`) is **not** about `Parameters`: **[verified]** `parameters.h` declares no `operator[]` in either `devel` or `branch_v1.8.x` (count 0). The NEWS line does not name the class, so its target is unconfirmed — but it cannot be `Parameters`.
- **`DofMap::dof_indices`** — no variadic/template change and no removal seen. **[verified]** the overload set (`elem`, `node`, `elem&`, `old_dof_indices`, `SCALAR_dof_indices`) is the same at `branch_v1.7.x` and `devel`; the range only *adds* `array_dof_indices(...)` and an extra `dof_indices(elem, ..., const std::vector<dof_id_type> &)` overload. `DofMap: add const DofConstraints accessor` and `Add a DofMap API to set need_full_sparsity_pattern` are additive.
- **`n_quadrature_points()`** — `1.8.0 -> 1.8.1`: "Move `n_quadrature_points()` to `FEAbstract`, fix it (#4137)". **[verified]** it is `virtual unsigned int n_quadrature_points() const = 0;` (pure) at `branch_v1.7.x` and non-pure with an implementation at `devel`. This *relaxes* the contract for FE subclasses (no longer required to implement it) and does not change call sites. `FEMContext` declares `n_quadrature_points` in neither version, so nothing to fix there.
- **`get_boundary_info()`** — no removal or signature change; it is the *recommended* replacement for the deprecated member. **[verified]** const and non-const overloads returning `BoundaryInfo &`. `BoundaryInfo::get_id_by_name(std::string_view)` and `MeshBase::get_id_by_name(std::string_view)` both still exist.
- **`get_id_by_name` name-lookup API itself** — unchanged; only its `regenerate_id_sets` input handling changed (see behavioural list).
- **FE `reinit` / `get_JxW` / `get_phi` / `get_dphi` / `get_xyz`** — no NEWS line in this range mentions these signatures. `FEAbstract::reinit(...)` is still `virtual` with the same shape. Nothing seen.
- **Quadrature / `QGauss`** — no NEWS mentions at all in this range; no rename, no added enum argument. (Only quadrature-adjacent lines: "Add nodal quadrature rules for Pyramids", "Add `clone` method to `QBase`", "Add `allow_nodal_pyramid_quadrature` flag" — all additive.)
- **`QoISet` / qoi indices** — **no NEWS lines in this range.** I found no mention, so I make no claim either way.
- **GetPot** — the only line in the whole range is `1.7.3 -> 1.8.0`: "Fix missing `<memory>` include in `getpot.h`", i.e. a build fix, not an API or include-path change.
- **`attach_assemble_function` / `attach_assemble_object` / `solve`** — no NEWS line changes these; **[verified]** `void attach_assemble_function (void fptr(EquationSystems &, const std::string &))` and `void attach_assemble_object (Assembly &)` are both present in `system.h`. `System::solve()`, `ImplicitSystem::solve() override`, `LinearImplicitSystem::solve() override` all present.
- **`System`/`EquationSystems` additive APIs** (safe, no action): "Add `System::has_constraint_object()`", "Add System matrix iterator APIs to match vector iterator APIs", "Add `CondensedEigenSystem::initialize_condensed_matrices()`, `have_condensed_dofs()`", "Allow belated `System::add_matrix()`", "Added `System::solve_for_unconstrained_dofs`", "`DofMap::process_mesh_constraint_rows`", "`DofMap::constrained_sparsity_construction()` option".
- **`NumericVector`/`SparseMatrix` additive APIs** (safe): `{NumericVector,SparseMatrix}::fuzzy_equals()`, `SparseMatrix::read_matlab()`, `col_start()/stop()`, `solver_package()`, `PetscMatrix::create_submatrix_nosort()`, `add_sparse_matrix()`, `PetscVector::pointwise_divide()`, `LumpedMassMatrix`, "Add matrixconvert app".
- **`MeshBase` additive APIs** (safe): `change_elemset_code()`, `get_elemset_codes()`, `active_unpartitioned_element_ptr_range()`, `reinit_ghosting_functors()`, `operator==`, `copy_constraint_rows()`, `n_constraint_rows()`, `print_constraint_rows()`, `elem_orders()`, `BoundaryInfo::operator==`, `BoundaryInfo::libmesh_assert_valid_multimaps()`, `MeshTools::Modification::orient_elements()`, `MeshTools::volume()`, `SimplexRefiner`, `Poly2TriTriangulator`, `STLIO`, `NonManifoldGhostingFunctor`/`SidesToElemMap`, `MeshfreeInterpolationFunction`.
- **Minimium library/compiler versions:** the only explicit requirement stated in this range is **C++17** ("C++17 support is now required", `1.7.3 -> 1.8.0`). On PETSc the NEWS is compatibility-progress, not a stated minimum: "Full compatibility with PETSc 3.17" (`1.7.0 -> 1.7.1`), "Full compatibility with PETSc 3.18 and 3.19" (`1.7.1 -> 1.7.3`), "PETSc 3.19 deprecated some convergence reasons", "Workaround for PETSc 3.18 losing %D support", "Convert `PETSC_RELEASE_LESS_THAN` macro usage to `PETSC_VERSION_LESS_THAN`". **No minimum PETSc version is stated**, and no ABI statement appears anywhere in the range.
- **Explicitly could not confirm (stated rather than guessed):**
  - "`ImplicitSystem::get_linear_solver()` now throws an error" — no throw found in the current implementation (see breaking list).
  - "Remove method deprecated in C++17" (`1.7.0 -> 1.7.1`) — the NEWS does not name the method.
  - "`inquire()` function signature change, move to anonymous namespace" (`1.6.2 -> 1.7.0`) — unnamed class/signature; reads internal.
  - "change in `InfFE::reinit(elem, pts):`" (`1.6.2 -> 1.7.0`) — the NEWS sentence is truncated after the colon, so the actual change is unstated.
  - "Macroify all the mesh iterators" (`1.7.3 -> 1.8.0`) — I did not verify whether public iterator *type names* changed; treat as worth a compile test if you name `MeshBase::element_iterator` etc. explicitly.
  - "Use `std::string_view` in place of `std::string` references" (`1.7.3 -> 1.8.0`) — this is why `get`/`have_parameter`/`get_id_by_name` now take `std::string_view` (**[verified]**). Ordinary calls with `std::string`/`const char *` are unaffected; only code taking member-function addresses or passing non-convertible types breaks.

---

## Quotes

All from `/home/wangchenyang/install/libmesh/NEWS` (line numbers on the left).

```
  2:   * Allow std::set Parameter instances.
  3:   * Avoid out of bounds access w/ operator[] and assert correctness.
  8:   * Make MeshBase::write() const.
 10:   * Don't look for 2D sides of a 2D element.
 19:   * Move n_quadrature_points() to FEAbstract, fix it (#4137)
 53:     - Add improved ASPECT_RATIO quality metric for Quad elements with unit tests
 55:     - Rename existing, generic ASPECT_RATIO metric to EDGE_LENGTH_RATIO
114:     - Add BoundaryInfo::renumber_id()
118:     - Allow child elements to have BoundaryInfo associated with them
124:     - Add option for removing global boundary ID in remove_id()
148:     - C++17 support is now required
162:     - Use std::make_unique everywhere, deprecate libmesh_make_unique
163:     - Use std::string_view in place of std::string references
164:     - Remove long-deprecated overload MeshTools::find_boundary_nodes()
165:     - Remove deprected UniquePtr workaround for std::unique_ptr
175:     - Require support for enum forward declarations
187:     - Change unsigned int to size_t in CouplingMatrix to prevent overflow for large number of variables (#3541)
205:     - Add a DofMap API to set need_full_sparsity_pattern
210:     - Fix missing <memory> include in getpot.h
221:     - PETSc 3.19 deprecated some convergence reasons
238:     - Pass unique_ptr by value to NumericVector::restore_subvector()
241:     - Add PetscMatrixBase<T>, base class for PetscMatrix<T> and PetscMatrixShellMatrix<T>
267:     - Deprecate RBParameters::get_parameters_map()
268:     - RBParameters: store std::vector<Real> instead of single Real
276:     - Replace NumericVector with std::vector in RBConstructionBase
294:     - Throw an error in EigenSystem::solve if we have constraints
297:     - Make System::is_initialized() const
300:     - Implement (Equation)System(s)::init() in terms of reinit_mesh()
308:     - Full compatibility with PETSc 3.18 and 3.19
312:     - Full compatibility with PETSc 3.17
313:     - Remove method deprecated in C++17
321:     - Act on id to name maps in BoundaryInfo::regenerate_id_sets
322:     - Selectively remove boundary information after mesh stitching
325:     - MeshBase reference member in BoundaryInfo changed to pointer
336:     - Add Elem::true_centroid(), Elem::vertex_average(), deprecate Elem::centroid()
348:     - Deprecate the SideEdge class
354:     - Deprecate the Side class
356:     - Stop building proxy sides by default
357:     - Set subdomain_id+p_level in all side/edge elems
362:     - Get exodusII.h include out of our headers
364:     - inquire() function signature change, move to anonymous namespace
373:     - Make ExodusII_IO_Helper::nodal_var_values a map
397:     - Use _active_vars in FEMContext calculations
398:     - Sort FEMContext::_active_vars
399:     - Optional active_vars in all FEMContext::FEMContext
400:     - Add active_vars member, accessor to FEMContext
415:     - change in InfFE::reinit(elem, pts):
439:     - Deprecate SparsityPattern::Build direct access
442:     - Allow belated System::add_matrix()
472:     - Include unpartitioned elements in MeshBase::subdomain_ids
478:     - Deprecate public MeshBase::boundary_info member
483:     - MeshFunction: use unique_ptr to manage PointLocator
494:     - System: store unique_ptrs in _vectors map
498:     - PointLocatorTree: use std::shared_ptr to manage tree
503:     - Officially deprecate UniquePtr
507:     - DofMap: add const DofConstraints accessor
528:     - DofMap::distribute_dofs() now returns the total number of DOFs across all procs
534:     - Re-enable check_dirichlet_bcid_consistency
535:     - Move linear_solver member from LinearImplicitSystem to ImplicitSystem
536:     - ImplicitSystem::get_linear_solver() now throws an error
537:     - UnsteadySolver: old_local_nonlinear_solution is now a std::shared_ptr
538:     - PeriodicBoundaries now stores unique_ptrs instead of dumb ptrs
543:     - ParsedFEMFunction: change "parsers" to a vector of unique_ptrs
571:     - Add virtual NumericVector::max_allowed_id() API.
573:     - NumericVector: Make indexing add/set operations toggle _is_closed
574:     - DistributedVector GHOSTED in serial can be SERIAL
603:     - Use SparsityPattern for n_nz in SparseMatrix
604:     - Make DofMap::build_sparsity() public
605:     - DofMap::get_sparsity_pattern() accessor
610:     - Remove SparsityPattern::Build::mesh
```

Version headings in the range read: `1.8.1 -> 1.8.2` (line 1), `1.8.0 -> 1.8.1` (15), `1.7.3 -> 1.8.0` (25), `1.7.1 -> 1.7.3` (306), `1.7.0 -> 1.7.1` (309), `1.6.2 -> 1.7.0` (315).

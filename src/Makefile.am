lib_LIBRARIES = libneko.a
libneko_a_SOURCES = \
	config/num_types.f90\
        common/structs.f90\
	common/parameters.f90\
	io/format/nmsh.f90\
	io/format/re2.f90\
	io/format/map.f90\
	io/format/rea.f90\
	common/log.f90\
	comm/comm.f90\
        mesh/curve.f90\
	comm/mpi_types.f90 \
	common/datadist.f90 \
	common/distdata.f90 \
	common/utils.f90\
	math/mxm_wrapper.F90\
	bcknd/cpu/mxm_std.f90\
	sem/speclib.f90 \
	math/math.f90\
	math/mathops.f90\
	math/fast3d.f90\
	math/tensor.f90\
	sem/space.f90\
	sem/dofmap.f90\
	sem/coef.f90\
	gs/gs_bcknd.f90\
	bcknd/cpu/gs_cpu.f90\
	bcknd/sx/gs_sx.f90\
	gs/gs_ops.f90\
	gs/gather_scatter.f90\
	mesh/entity.f90\
	mesh/point.f90\
	mesh/element.f90\
	math/ax.f90\
	ksp/precon.f90\
	ksp/krylov.f90\
	ksp/cg.f90\
        ksp/bicgstab.f90\
	mesh/quad.f90\
	mesh/hex.f90\
	mesh/tet.f90\
	mesh/tri.f90\
	adt/htable.f90\
	adt/uset.f90\
	adt/stack.f90\
	adt/tuple.f90\
	mesh/zone.f90\
	mesh/mesh.f90\
	field/mesh_field.f90\
	mesh/tet_mesh.f90\
	field/field.f90\
	bcknd/sx/sx_cdtp.f90\
	bcknd/sx/sx_conv1.f90\
	bcknd/sx/sx_dudxyz.f90\
	bcknd/sx/sx_opgrad.f90\
	math/operators.f90\
	bcknd/cpu/opr_cpu.f90\
	bcknd/sx/opr_sx.f90\
	bcknd/xsmm/opr_xsmm.F90\
	io/generic_file.f90\
	io/map_file.f90\
	io/re2_file.f90\
	io/rea_file.f90\
	io/fld_file.f90\
	io/vtk_file.f90\
	io/nmsh_file.f90\
	io/file.f90\
	io/output.f90\
	io/fluid_output.f90\
	io/sampler.f90\
	bc/bc.f90\
	bc/dirichlet.f90\
	bc/wall.f90\
	bc/inflow.f90\
	bc/usr_inflow.f90\
	bc/facet_normal.f90\
	bc/symmetry.f90\
	ksp/gmres.f90\
	ksp/pc_jacobi.f90\
	ksp/pc_identity.f90\
	bcknd/sx/pc_jacobi_sx.f90\
	common/source.f90\
	common/abbdf.f90\
	config/neko_config.f90\
	case.f90\
	common/user_intf.f90\
	solver/fluid_method.f90\
	solver/fluid_plan4.f90\
	solver/fluid_plan1.f90\
	solver/fluid_schemes.f90\
	simulation.f90\
	math/ax_helm.f90\
	bcknd/sx/ax_helm_sx.f90\
	bcknd/xsmm/ax_helm_xsmm.F90\
	solver/projection.f90\
	comm/parmetis.F90\
	comm/redist.f90\
	ksp/pc_hsmg.f90\
	math/fdm.f90\
	math/schwarz.f90\
	neko.f90

if ENABLE_PARMETIS
libneko_a_SOURCES += comm/parmetis_wrapper.c
endif

bin_PROGRAMS = neko
neko_SOURCES = $(libneko_a_SOURCES) driver.f90

pkginclude_HEADERS =  $(wildcard *.mod)

CLEANFILES =*.mod

if ENABLE_MAKEDEPF90
depend:
	$(MAKEDEPF90) $(libneko_a_SOURCES) driver.f90  > .depends
endif

include	.depends

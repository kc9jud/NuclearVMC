$(eval $(begin-module))

################################################################
# unit definitions
################################################################

module_units_h := units random_buffer av18pot
#module_units_cpp-h := wavefunction
module_units_f := av18pot
module_programs_cpp := random_buffer_test
# module_programs_f :=
# module_generated :=

################################################################
# library creation flag
################################################################

$(eval $(library))

$(eval $(end-module))

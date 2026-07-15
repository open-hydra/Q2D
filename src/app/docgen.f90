program Q2D_docgen
    ! Q2D docgen driver: builds the registry of all INI parameters by invoking
    ! every Register_* routine, then generates `docs/user/registry.md` via
    ! `reg%generate_markdown`.
    !
    use Q2D_Read_Sim_Param, only: Register_Sim_Param
    use Q2D_Read_IO,        only: Register_IO_Fields, Register_Probes, Register_LineProbes
    use Q2D_Read_Numerics,  only: Register_Numerics
    use Q2D_Read_Physics,   only: Register_Physics
    use Q2D_Read_Sources,   only: Register_Sources
    use Q2D_Registry_m,     only: reg
    use Q2D_Config_Types_m, only: obj_chemistry
    implicit none

    ! Placeholder allocation of the chemistry tolerance arrays (stub nsc=3 +
    ! 1 slot for temperature): required because Register_Physics writes into the
    ! obj_chemistry%RT/AT sections.
    allocate(obj_chemistry%RT(1:3+1), obj_chemistry%AT(1:3+1))

    ! Build the registry entries from all six Q2D registration sources.
    call Register_Sim_Param()
    call Register_IO_Fields()
    call Register_Probes(1, 'probe-placeholder')  ! Register one probe with a placeholder name
    call Register_LineProbes(1, ['lineprobe-placeholder']) 
    call Register_Numerics(2)                      ! placeholder nmgl=2 (multigrid levels)
    call Register_Physics()
    call Register_Sources()                        ! geometric source-term parameters

    call reg%generate_markdown('docs/user/registry.md')

end program Q2D_docgen

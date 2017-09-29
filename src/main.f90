#ifdef _WIN32

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#DEFINE DASH "\"
#DEFINE COPY "copy"
#DEFINE MD "md"
#DEFINE RMDIR "RMDIR /S /Q"
#DEFINE RM "del"
#DEFINE RENAME "MOVE /Y"
#DEFINE SH "BAT"
#DEFINE EXE ".exe"
#DEFINE NULL " >NUL"


#else

#define STRINGIFY(x)#x
#define TOSTRING(x) STRINGIFY(x)

#DEFINE DASH "/"
#DEFINE COPY "cp"
#DEFINE MD "mkdir"
#DEFINE RMDIR "rm -r"
#DEFINE RM "rm"
#DEFINE RENAME "mv"
#DEFINE SH "sh"
#DEFINE EXE ""
#DEFINE NULL ""


#endif
!######################################################################

program AlphaImpute
    ! The program is intended to be run numerous times
    ! The first time AlphaImpute is run in order to calculate Genotype probabilities (GenoProb)
    ! The second time, it is run to phase the individuals
    ! Finally, AlphaImpute should be run to impute genotypes.
    ! RestartOptions handle this situation, so:
    !   * inputParams%restartOption=0 => Passes through the whole process: GenoProb, Phasing and Imputing
    !   * inputParams%restartOption=1 => Makes only GenoProbReadInParameterFile
    !   * inputParams%restartOption=2 => Makes GenoProb and Phasing
    !   * inputParams%restartOption=3 => Makes only Imputation. This implies AlphaImpute has to be run
    !                           already in order to get GenoProb done or Genotype Probabilities
    !                           Genotype Probabilities have to be edited by hand
    !   * inputParams%restartOption=4 =>
    use Global
    use AlphaImputeModule
    use AlphaImputeSpecFileModule
    use alphaFullChromModule
    ! use alphaFullChromModule

    implicit none

    character(len=4096) :: cmd, SpecFile

    if (Command_Argument_Count() > 0) then
        call get_command_argument(1,cmd)
        if (cmd(1:2) .eq. "-v") then
            call PrintVersion
            call exit(0)
        end if
    end if

    if (Command_Argument_Count() > 0) then
        call Get_Command_Argument(1,SpecFile)
    else
        specfile="AlphaImputeSpec.txt"
    end if



    call Titles


    ! use default input, TODO this can be changed
    allocate(defaultInput)
    call defaultInput%ReadInParameterFile(SpecFile)


    if (defaultInput%plinkinputfile == "") then
        call runPlink(defaultInput%plinkinputfile, defaultInput, runAlphaImpute)

    else
        call runAlphaImpute(defaultInput)
    endif


 

end program AlphaImpute


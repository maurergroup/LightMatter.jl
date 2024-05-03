using DelimitedFiles,Unitful

abstract type Method end
abstract type ETTM <: Method end
abstract type TTM <: Method end

@kwdef struct simulationsettings{M} <:ETTM{M}

end

@kwdef struct simulationsettings <:TTM

end

function Units(dict::Dict)
    dict["g"]=ustrip(uconvert(u"eV/nm^3/fs/K",dict["g"]u"W/m^3/K"))
    dict["gamma"]=ustrip(uconvert(u"eV/nm^3/K^2",dict["gamma"]u"J/m^3/K^2"))
    dict["ExtCof"]=ustrip(uconvert(u"nm",dict["ExtCof"]u"m"))
    dict["AtomDens"]=ustrip(uconvert(u"nm^-3",dict["AtomDens"]u"m^-3"))
    dict["FWHM"]=ustrip(uconvert(u"fs",dict["FWHM"]u"s"))
    dict["Fluence"]=ustrip(uconvert(u"eV/nm^2",dict["Fluence"]u"J/m^2"))
    dict["LaserOff"]=ustrip(uconvert(u"fs",dict["LaserOff"]u"s"))
    dict["RTKappa"]=ustrip(uconvert(u"eV/nm/K/fs",dict["RTKappa"]u"W/m/K"))
    dict["SimEnd"]=ustrip(uconvert(u"fs",dict["SimEnd"]u"s"))
    dict["Length"]=ustrip(uconvert(u"nm",dict["Length"]u"m"))
    dict["dz"]=ustrip(uconvert(u"nm",dict["dz"]u"m"))
    dict["Plasma"]=dict["Plasma"]/(6.582e-1*2*pi) #Plasma/2pihbar
end

function generate_inputs(file::String)
    IO = readdlm(file,'=',comments=true,comment_char='#')
    for (x,y) in enumerate(IO[:,2])
        if typeof(y) == SubString{String}
            IO[x,2]=strip(IO[x,2],' ')
        end
    end
    Input=Dict(IO[i,1]=>IO[i,2] for i in 1:size(IO,1))
    if Input["Type"]=="E2TM"
        Units(Input)
    end
    return Input
end



function define_simulation_settings()

end


function define_simulation_settings(dict::Dict)

end

function Unit_conversion(dict::Dict)
    dict["g"]=ustrip(uconvert(u"eV/nm^3/fs/K",dict["g"]))
    dict["Gamma"]=ustrip(uconvert(u"eV/nm^3/K^2",dict["Gamma"]))
    dict["ExtCof"]=ustrip(uconvert(u"nm",dict["ExtCof"]))
    dict["AtomDens"]=ustrip(uconvert(u"nm^-3",dict["AtomDens"]))
    dict["FWHM"]=ustrip(uconvert(u"fs",dict["FWHM"]))
    dict["Fluence"]=ustrip(uconvert(u"eV/nm^2",dict["Fluence"]))
    dict["LaserOff"]=ustrip(uconvert(u"fs",dict["LaserOff"]))
    dict["RTKappa"]=ustrip(uconvert(u"eV/nm/K/fs",dict["RTKappa"]))
    dict["SimEnd"]=ustrip(uconvert(u"fs",dict["SimEnd"]))
    dict["Length"]=ustrip(uconvert(u"nm",dict["Length"]))
    dict["dz"]=ustrip(uconvert(u"nm",dict["dz"]))
    dict["Plasma"]=ustrip(dict["Plasma"])
end

function generate_inputs(file::String)
    IO = readdlm(file,':',comments=true,comment_char='#')
    for (x,y) in enumerate(IO[:,2])
        if typeof(y) == SubString{String}
            IO[x,2]=strip(IO[x,2],' ')
        end
        if IO[x,3] != ""
            IO[x,2] = y* uparse(IO[x,3])
        end
    end
    Input=Dict(IO[i,1]=>IO[i,2] for i in 1:size(IO,1))
    Unit_conversion(Input)
    println(Input["Plasma"])
end
#Test=generate_inputs("InputFiles/Au_Input.txt")
generate_inputs("InputFiles/Au_Input.txt")
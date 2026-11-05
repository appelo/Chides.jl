using Pkg
Pkg.activate(@__DIR__)
cd(@__DIR__) do
    Pkg.develop(path="../..")
end
Pkg.instantiate()
println("Chapter 2 environment ready. Build the Fortran library with make -C src/chides_fortran.")

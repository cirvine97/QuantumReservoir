FROM julia:1.10.4

WORKDIR /app

COPY . /app

RUN julia -e 'using Pkg; Pkg.add(["LinearAlgebra", "QuantumInformation", "Plots", "ProgressBars", "DataFrames", "CSV", "Random", "JSON"])'

CMD ["julia", "generate_experimental_data.jl"]

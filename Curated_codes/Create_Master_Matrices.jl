######################################################################## # Codes for the paper:
# An Anthropocene Map of Genetic Diversity
#
# Andreia Miraldo, Sen Li, Michael K. Borregaard, Alexander Floréz-Rodriguéz, Shyam Gopalakrishnan, Mirneza Risvanovic, Zhiheng Wang, Carsten Rahbek, Katharine A. Marske & David Nogués-Bravo
#
# Submitted to Science, 2016
# Code in this file by Sen Li and Michael K. Borregaard
#########################################################################

# Load the necessary library and functions
include("Compare_Pairwise_Function.jl")
using DataFrames

"""
A function to create master data matrices that are used to compute genetic diversity, assess data quality and do sensitivity analyses.

**Parameters**
* 'foldername' : The name of of a folder containing the data files. There must be files of two types (file extension 'fasta' and file extension 'coords') with  the same filename, e.g. the species names (e.g. folder contents could be 'Bufo_bufo.fasta, Bufo_bufo.coords, Rana_arvalis.fasta, Rana_arvalis.coords', etc.). The .fasta files contain the sequences as an m x n integer matrix, where m is the number of sequences and n is the length of the alignments. The .coords files contain the geographic coordinates of the sequences, as an m x 2 floating point matrix with latitude in the first column and longitude in the second.  TODO: A folder with the anthrome and biome of each sequence
"""
function create_master_matrices(foldername::ASCIIString)
    species_list = [x[1:(end-6)] for x in filter(st->contains(st, ".fasta"), readdir(foldername))]      #identify unique file names ignoring the extension
    num_files = size(species_list,1)

    equalarea = latbands = biomes = anthrome = gridcells = DataFrame(species = UTF8String[], cell = UTF8String[], seq1 = Int[], seq2 = Int[], length_seq1 = Int[], length_seq2 = Int[], overlap = Float64[], commons = Int[], num_per_bp = Float64[])   # Pre-initialize the DataFrame to ensure correct element types

    for iter_file = 1:num_files
        # read in the sequence data for one species
        file_name_seq = joinpath(foldername, species_list[iter_file] * ".fasta")
        species_seqs = readdlm(file_name_seq, Int)
        length_seq = size(species_seqs,2)

        # read the coordinates for one species
        file_name_coords = joinpath(foldername, species_list[iter_file] * ".coords")
        species_coords = readdlm(file_name_coords)

        #GRID CELLS
        #apply a 4x4 lat-long degree grid to the coordinates
        coords_grid = hcat(ceil(Int, (species_coords[:,1] + 2) ./ 4) .* 4 .- 4, ceil(Int, species_coords[:,2] ./ 4) .* 4 .- 2)

        #calculate the summary statistic for each grid cell
        sitenames = mapslices(x -> "$(x[1])_$(x[2])", coords_grid, 2)
        gridcells = vcat(gridcells, calcspecies(species_list[iter_file], species_seqs, sitenames))

        #ANTHROMES
        file_name_anthrome = joinpath("biome_anthrome_latband", foldername, foldername[end-3:end] * "_anthrome_matlab",species_list[iter_file] * ".coords")
        anthro = floor(Int, readdlm(file_name_anthrome)[:,3])
        # calculate the summary statistic for each grid cell
        anthromes = vcat(anthromes, calcspecies(species_list[iter_file], species_seqs, anthro))

        #BIOMES
        file_name_biome = joinpath("biome_anthrome_latband", foldername, foldername[end-3:end] * "_biome_matlab",species_list[iter_file] * ".coords")
        bio = floor(Int, readdlm(file_name_biome)[:,3])
        biomes = vcat(biomes, calcspecies(species_list[iter_file], species_seqs, bio))

        #LATITUDINAL BANDS
        latband = 10 * floor(Int,species_coords[:,1]/10)
        latbands = vcat(latbands, calcspecies(species_list[iter_file], species_seqs, latband))

        #EQUAL AREA GRID CELLS
        file_name_equal = joinpath("biome_anthrome_latband", foldername, foldername[end-3:end] * "_equalarea_matlab",species_list[iter_file] * ".coords")
        equ = floor(Int, readdlm(file_name_equal)[:,3])
        equ = [UTF8String("$a") for a in equ]
        equalarea = vcat(equalarea, calcspecies(UTF8String(species_list[iter_file]), species_seqs, equ))

        println(iter_file)

    end
    writetable("pairwise_gridcells_$(foldername).csv", gridcells)
    writetable("pairwise_anthromes_$(foldername).csv", anthromes)
    writetable("pairwise_biomes_$(foldername).csv", biomes)
    writetable("pairwise_latbands_$(foldername).csv", latbands)
    writetable("pairwise_equalarea_$(foldername).csv", equalarea)
end

"""
A function to calculate the summary statistic for all sites (e.g. grid cell or biome) for one species.

**Parameters**
* 'species':        A string with the name of the species
* 'species_seqs':   A matrix where the rows are aligned genetic sequences, and columns are loci. Basepairs must be coded as 1, 2, 3 or 4, or with a 0 signifying that the locus is absent from the alignment.
* 'sitenames':      A vector of strings with the names of the sites
"""
function calcspecies(species::UTF8String, species_seqs::Matrix{Int}, sitenames::Vector{UTF8String})
    res = DataFrame(species = UTF8String[], cell = UTF8String[], seq1 = Int[], seq2 = Int[], length_seq1 = Int[], length_seq2 = Int[], overlap = Float64[], commons = Int[], num_per_bp = Float64[])
    uniq_grid = unique(sitenames, 1)
    for iter_grid = 1:size(uniq_grid,1) # Go through each site in turn
        seqs_grid = species_seqs[findin(sitenames, uniq_grid[iter_grid, :]), :]

        tot_muts = compare_pairwise(seqs_grid)
        lns = size(tot_muts, 1)
        tmp = DataFrame(species = rep(species, lns), cell = rep(vec(uniq_grid[iter_grid, :]), lns)) #Expand species and cell names to the length of the resulting DataFrame
        tot_muts = hcat(tmp, tot_muts)

        res = vcat(res, tot_muts)
    end

    for sym in [:seq1, :seq2, :length_seq1, :length_seq2, :overlap, :commons, :num_per_bp]
        res[sym][res[sym] .< 0.] = NA   #Replace empty values with an NA identifier
    end
    res
end

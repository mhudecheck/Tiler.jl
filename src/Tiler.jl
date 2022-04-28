module Tiler
    using Sentinel
    using ArchGDAL
    using TiledIteration
    using Images
    using ColorSchemes
    using CUDA
    using Base.Threads

    export resizeTile, saveTiles, mapTile, formatTIF, mapTIF, loadTIF, warpTIF

    function loadTIF(tif; band=1)
        @show "loadTIF()"
        ArchGDAL.read(tif) do dataset
            @show "Reading TIF with ArchGDAL"
            resBand = ArchGDAL.getband(dataset, band)
            b = ArchGDAL.read(resBand)
            @show "TIF Loaded loadTIF()"
            #return b
            return rotl90(reverse(b, dims=2))
        end
    end

    #=function formatTIF(tif::String; axis::Real=0, rotate::Bool=false)
        if axis == 0
            gdalString = ["-te_srs", "EPSG:4326", "-te", "-180", "-85.051129", "180", "85.051129", "-t_srs", "EPSG:3857"]
        else
            gdalString = ["-te_srs", "EPSG:4326", "-te", "-180", "-85.051129", "180", "85.051129", "-t_srs", "EPSG:3857", "-ts", string(axis[2]), string(axis[1])]
        end
        ArchGDAL.read(tif) do test
            ArchGDAL.gdalwarp([test], gdalString) do dataset
                numBands = ArchGDAL.nraster(dataset)
                origBand = ArchGDAL.getband(dataset, 1)
                b1 = ArchGDAL.read(origBand)
                if rotate != false
                    b1 = rotl90(reverse(b1,dims=2))
                end
                if numBands != 1
                    varType = ArchGDAL.pixeltype(origBand)
                    axisWidth = ArchGDAL.width(origBand)
                    axisHeight = ArchGDAL.height(origBand)
                    if rotate != false
                        outputArray = Array{Float32}(undef, axisHeight, axisWidth, numBands)
                    else
                        outputArray = Array{Float32}(undef, axisWidth, axisHeight, numBands)
                    end
                    outputArray[:,:,1] = b1 ./ 255
                    for i in 2:numBands
                        band = ArchGDAL.getband(dataset,i)
                        b2 = ArchGDAL.read(band)
                        if rotate != false
                            b2 = rotl90(reverse(b2,dims=2))
                        end
                        outputArray[:,:,i] = b2 ./ 255
                    end
                    return outputArray
                else
                    return b1
                end
            end
        end
    end=#

    function warpTIF(tifName::String; sourceDirectory::String="./", tempDirectory::String="./", numThreads::Real=1)
        fullTif = sourceDirectory * "/" * tifName
        fullOutputTif = tempDirectory * "/" * "warped_" * tifName
        if isfile(fullOutputTif)
            rm(fullOutputTif)
        end
        #if isfile(fullOutputTif) != true
        aggr = `gdalwarp  -te_srs EPSG:4326 -te -180 -85.051129 180 85.051129 -t_srs EPSG:3857 -co NUM_THREADS=$numThreads -wo NUM_THREADS=$numThreads $fullTif $fullOutputTif`
        run(aggr)
        #end
    end

    function formatTIF(tifName::String; sourceDirectory::String="./", tempDirectory::String="./", axis = nothing, transpose = false, numThreads = 1, removeTif = true)
        warpTIF(tifName; sourceDirectory=sourceDirectory, tempDirectory=tempDirectory, numThreads=numThreads)
        fullOutputTif = tempDirectory * "/" * "warped_" * tifName
        rasterIn = ArchGDAL.readraster(fullOutputTif)
        rType = typeof(rasterIn[1,1,1])
        rWidth = size(rasterIn, 1)
        rHeight = size(rasterIn, 2)
        nBands = size(rasterIn, 3)
        #println(nBands)
        #println(rWidth)
        #println(rHeight)
        returnArray = Array{Float16}(undef, nBands, rHeight, rWidth)
        for i in 1:nBands
            inArray = rasterIn[:, :, i]
            testArray = convert(Array{Float16}, inArray)
            #finalArray = cudaCleanTranspose(testArray, 32)
            returnArray[i, :, :] = testArray'
        end
        if removeTif == true
            rm(fullOutputTif)
        end
        return returnArray
    end

    #=function mapTIF(tif::String, outputName::String; outputDirectory::String=".", scheme::String="inferno", min::Real=0, max::Real=1, type=UInt8)
        @show scheme
        ArchGDAL.read(tif) do dataset
            # Read Information
            band = ArchGDAL.getband(dataset, 1)
            data = ArchGDAL.read(band)
            bandType = ArchGDAL.pixeltype(band)
            bandWidth = ArchGDAL.width(dataset)
            bandHeight = ArchGDAL.height(dataset)
            ref = ArchGDAL.getproj(dataset)
            geotransform = ArchGDAL.getgeotransform(dataset)
    
            # Apply Color Map
            convertedData = get(colorschemes[Symbol(scheme)], data, (convert(bandType, min), convert(bandType, max)))
            chanView = channelview(convertedData)
            finDat = permuteddimsview(chanView, (2,3,1))
    
            # Create New Raster
            mkpath(outputDirectory)
            output = outputDirectory * "/" * outputName
            ArchGDAL.create(output, driver=ArchGDAL.getdriver("GTiff"), width=bandWidth, height=bandHeight, nbands=3, dtype=type) do raster
                ArchGDAL.setgeotransform!(raster, geotransform)
                ArchGDAL.setproj!(raster, ref)
                ArchGDAL.write!(raster, finDat[:,:,1]*255, 1)
                ArchGDAL.write!(raster, finDat[:,:,2]*255, 2)
                ArchGDAL.write!(raster, finDat[:,:,3]*255, 3)
            end
            return convertedData
        end
    end=#

    # Applies Color Scheme to GeoTIFF and Saves Result
    function mapTIF(tif::String, outputName::String; outputDirectory::String=".", scheme::String="inferno", min::Real=0, max::Real=1, type=UInt8, scale::Real=1)
        ArchGDAL.read(tif) do dataset
            # Read Information
            band = ArchGDAL.getband(dataset, 1)
            data = ArchGDAL.read(band)
            bandType = ArchGDAL.pixeltype(band)
            bandWidth = ArchGDAL.width(dataset)
            bandHeight = ArchGDAL.height(dataset)
            ref = ArchGDAL.getproj(dataset)
            geotransform = ArchGDAL.getgeotransform(dataset)

            # Apply Color Map
            convertedData = get(colorschemes[Symbol(scheme)], data, (convert(bandType, min), convert(bandType, max)))
            chanView = channelview(convertedData)
            finDat = permuteddimsview(chanView, (2,3,1))

            # Create New Raster
            mkpath(outputDirectory)
            output = outputDirectory * "/" * outputName
            ArchGDAL.create(output, driver=ArchGDAL.getdriver("GTiff"), width=bandWidth, height=bandHeight, nbands=3, dtype=type) do raster
                ArchGDAL.setgeotransform!(raster, geotransform)
                ArchGDAL.setproj!(raster, ref)
                if scale != 1
                    ArchGDAL.write!(raster, finDat[:,:,1]*scale, 1)
                    ArchGDAL.write!(raster, finDat[:,:,2]*scale, 2)
                    ArchGDAL.write!(raster, finDat[:,:,3]*scale, 3)
                else
                    ArchGDAL.write!(raster, finDat[:,:,1], 1)
                    ArchGDAL.write!(raster, finDat[:,:,2], 2)
                    ArchGDAL.write!(raster, finDat[:,:,3], 3)
                end
            end
            return convertedData
        end
    end

    function mapTile(tif::String; sourceDirectory::String="./", gpu::Bool=false, tile::Real=256, min::Real=0, device=0, max=8, filter=false, alpha = 1.0, kernel = Kernel.gaussian(3), outputDirectory=".", scale = false, scaleLevel = 0, formatted=false, numThreads::Real=1, smoothing::String="linear", type=UInt8)
        mkpath(outputDirectory)
        if formatted == false
            img = formatTIF(tif; sourceDirectory=sourceDirectory, transpose=true, removeTif=false, numThreads = numThreads)
        else
            img = loadTIF(tif)
        end
        type = 0
        if size(img, 3) > 2 # RGB o
            dims = 4

            img = RGBA{N0f8}.(colorview(RGB, img), alpha)
            if filter == true
                img = RGBA{N0f8}.(imfilter(img, kernel))
            end
            type = eltype(img)
            imgTuple = reinterpret(NTuple{dims, UInt8}, img)
        elseif scale == true # Grayscale & Adjust Vals
            imgScaled = img .* scaleLevel
            imgTuple = convert(typeof(img), imgScaled)
        else # Grayscale
            for i in min:max
                #@time resizeTile(img, i; axis=tile, saveType=type, gpu=gpu, outputDirectory=outputDirectory)
            end
        end
        #for i in min:max
            #@info "Zoom level: $i"
        @time resizeTile(imgTuple, min, max; axis=tile, device=device, saveType=type, gpu=gpu, outputDirectory=outputDirectory, smoothing=smoothing)
        #end
        GC.gc()
        return 
    end

    function resizeTile(raster::AbstractArray, zoomMin::Integer=0, zoomMax::Integer=8; gpu::Bool=true, axis::Integer=256, saveType=0, outputDirectory::String=".", interpolation::Bool=true, smoothing::String="linear", device=0)
        conv = 0
        if gpu == true & CUDA.functional() == true
            if isa(raster, CuArray) == false
                raster = CuArray(raster)
                conv = 1
            end
        end
        for i in zoomMin:zoomMax
            @info "Zoom level: $i"
            targetDim = axis * 2^i
            @time outputRaster = resizeRaster(raster, targetDim, targetDim; gpu=gpu, device=device, interpolation=interpolation, smoothing=smoothing);
            @time saveTiles(outputRaster, axis; zoom=i, type=saveType, outputDirectory=outputDirectory)
            outputRaster = nothing
        end
        if conv == 1
            @time raster = Array(raster)
        end
    end

    function saveTiles(img::AbstractArray, axis::Real; zoom=0, type=0, outputDirectory::String=".", xCoord::Real=1, yCoord::Real=1, yMax::Real=1)
        tileinds_all = collect(TileIterator(axes(img), (axis, axis)))
        @info tileinds_all, 1:length(tileinds_all), axes(img), axis
        #k = sqrt(length(tileinds_all))

        Threads.@threads for j = 1:length(tileinds_all)
            y = (ceil(Int, (tileinds_all[j][1][2] / axis)) - 1) + (2^(zoom) * yCoord) - 2^(zoom)
            x = (ceil(Int, (tileinds_all[j][2][2] / axis)) - 1) + (2^(zoom) * xCoord) - 2^(zoom)
            tileAxis = tileinds_all[j]
            tile = img[tileAxis...]
            @info "Saving tiles: $zoom, $x, $y, $xCoord, $yCoord"
            if type != 0
                tile = reinterpret(type, tile)
                save("$outputDirectory/$zoom/$x/$y.png", tile)
            else
                tile = reinterpret(N0f16, tile) # Converts grayscale to normed UInt16
                save("$outputDirectory/$zoom/$x/$y.png", tile)
            end
        end
    end
end





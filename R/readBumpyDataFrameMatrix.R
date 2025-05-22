#' Read a BumpyDataFrameMatrix from disk
#'
#' Read a \linkS4class{BumpyDataFrameMatrix} from its on-disk representation.
#'
#' @param path String containing a path to a directory, itself created using the \code{\link{saveObject}} method for \linkS4class{BumpyDataFrameMatrix} objects.
#' @param metadata Named list of metadata for this object, see \code{\link{readObjectFile}} for details.
#' @param ... Further arguments passed to internal \code{\link{altReadObject}} calls.
#'
#' @return A \linkS4class{BumpyDataFrameMatrix} object.
#'
#' @author Aaron Lun
#' @examples
#' # Mocking up a BumpyMatrix.
#' library(BumpyMatrix)
#' library(S4Vectors)
#' df <- DataFrame(x=runif(100), y=runif(100))
#' f <- factor(sample(letters[1:20], nrow(df), replace=TRUE), letters[1:20])
#' out <- S4Vectors::split(df, f)
#' mat <- BumpyMatrix(out, c(5, 4))
#'
#' # Saving it:
#' tmp <- tempfile()
#' saveObject(mat, tmp)
#' 
#' # Reading it:
#' readBumpyDataFrameMatrix(tmp)
#'
#' @export
#' @aliases loadBumpyDataFrameMatrix
#' @importFrom BumpyMatrix BumpyMatrix
#' @importFrom IRanges PartitioningByWidth
#' @importFrom Matrix sparseMatrix
#' @importFrom BiocGenerics relist
#' @import alabaster.base
readBumpyDataFrameMatrix <- function(path, metadata, ...) {
    concatenated <- altReadObject(file.path(path, "concatenated"), ...)

    fhandle <- H5Fopen(file.path(path, "partitions.h5"), "H5F_ACC_RDONLY")
    on.exit(H5Fclose(fhandle), add=TRUE, after=FALSE)
    ghandle <- H5Gopen(fhandle, "bumpy_data_frame_array")
    on.exit(H5Gclose(ghandle), add=TRUE, after=FALSE)

    lengths <- h5_read_vector(ghandle, "lengths")
    dfl <- relist(concatenated, PartitioningByWidth(x = lengths))

    dimnames <- list(NULL, NULL)
    if (h5_object_exists(ghandle, "names")) {
        nhandle <- H5Gopen(ghandle, "names")
        on.exit(H5Gclose(nhandle), add=TRUE, after=FALSE)
        for (i in seq_along(dimnames)) {
            nm <- as.character(i - 1L)
            if (h5_object_exists(nhandle, nm)) {
                dimnames[[i]] <- h5_read_vector(nhandle, nm)
            }
        }
    }

    dims <- h5_read_vector(ghandle, "dimensions")
    if (dims[1] * dims[2] > length(dfl)) {
        ihandle <- H5Gopen(ghandle, "indices")
        rows <- h5_read_vector(ihandle, "0")
        cols <- h5_read_vector(ihandle, "1")
        mat <- sparseMatrix(i=rows + 1L, j=cols + 1L, x=seq_along(dfl), dims=dims, dimnames=dimnames)
    } else {
        mat <- matrix(seq_along(dfl), dims[1], dims[2], dimnames=dimnames)
    }

    BumpyMatrix(x = dfl, proxy = mat)
}

##########################
##### OLD STUFF HERE #####
##########################

#' @export
loadBumpyDataFrameMatrix <- function(assay.info, project) { 
    concat.meta <- acquireMetadata(project, assay.info$bumpy_data_frame_matrix$concatenated$resource$path)
    concatenated <- .loadObject(concat.meta, project)

    dimnames <- list(NULL, NULL)
    dimnames_names <- c("row_names", "column_names")
    for (d in 1:2) {
        if (!is.null(dim.info <- assay.info$bumpy_matrix[[dimnames_names[d]]])) {
           dim.meta <- acquireMetadata(project, dim.info$resource$path)
            dimnames[[d]] <- .loadObject(dim.meta, project)[,1]
        }
    }

    group.path <- acquireFile(project, assay.info$path)
    groupings <- .quickReadCsv(group.path, 
        expected.columns=c(row="integer", column="integer", number="integer"),
        expected.nrows=assay.info$bumpy_matrix$object_count,
        row.names=isTRUE(assay.info$bumpy_matrix$object_names),
        compression=assay.info$bumpy_matrix$compression
    )

    dfl <- relist(concatenated, PartitioningByWidth(x = groupings$number))
    if (!is.null(rownames(groupings))) {
        names(dfl) <- rownames(groupings)
    }

    dims <- assay.info$array$dimensions
    if (dims[[1]] * dims[[2]] > nrow(groupings)) {
        mat <- sparseMatrix(i=groupings$row, j=groupings$column, x=seq_len(nrow(groupings)), dims=as.integer(dims), dimnames=dimnames)
    } else {
        mat <- matrix(0L, dims[[1]], dims[[2]], dimnames=dimnames)
        mat[as.matrix(groupings[,c("row", "column")])] <- seq_len(nrow(groupings))
    }

    BumpyMatrix(x = dfl, proxy = mat)
}

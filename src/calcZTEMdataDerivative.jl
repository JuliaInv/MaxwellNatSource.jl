
export calcZTEMdata, calcFullZTEMDerivative

function calcZTEMdata( DD::Array{Complex128,2} )  # (ndata,2) E and H data from polarizations 1 and 2
    # Data must be ordered [ hx,hy,hz, ... ]
    # Calculate: 

    #          ( -Hy1*Hz2 + Hy2*Hz1 )
    # ( Tx )   (  Hx1*Hz2 - Hx2*Hz1 )
    # ( Ty ) = ----------------------
    #             Hx1*Hy2 - Hx2*Hy1

    fieldsPerData = 3

    ndata = size(DD,1)   
    nr = div(ndata, fieldsPerData)

    zdpred = Array{Complex128}(2,nr)

    id = 0

    for ir = 1:fieldsPerData:ndata
        Hx1 = DD[ir,   1]
        Hy1 = DD[ir+1, 1]
        Hz1 = DD[ir+2, 1]

        Hx2 = DD[ir,   2]
        Hy2 = DD[ir+1, 2]
        Hz2 = DD[ir+2, 2]

        HH = Hx1*Hy2 - Hx2*Hy1

        id += 1
        zdpred[1,id] = ( -Hy1*Hz2 + Hy2*Hz1 ) / HH
        zdpred[2,id] = (  Hx1*Hz2 - Hx2*Hz1 ) / HH
    end 

    return zdpred
end

#--------------------------------------------------------------------------

function calcFullZTEMDerivative( DD::Array{Complex128,2} )  # (nr*4,2) E and H data from polarizations 1 and 2
    # Data must be ordered [ hx,hy,hz, ... ]
    # Calculate: T

    # Output is a 2-block diagonal matrix.
    # In the derivative matrix, for each receiver, the 4 rows are ordered:
    #    [  Txr Txi  Tyr Tyi ]
    # and the 2*6 columns are ordered:
    #    [ Hx1r Hx1i  Hy1r Hy1i  Hz1r Hz1i ]
    #    [ Hx2r Hx2i  Hy2r Hy2i  Hz2r Hz2i ]

    fieldsPerData = 3

    ndata = size(DD,1) 
    if mod(ndata,fieldsPerData) != 0
        error("# of data should be divisible by 3.")
    end

    nr = div(ndata, fieldsPerData)  # number of data locations

    nrow =  4
    ncol =  6 
    dZdD = spzeros(Float64, nr*nrow, nr*ncol*2)

    j = 1
    for ir = 1:fieldsPerData:ndata
        deriv1,deriv2 = ZTEMderivs( DD[ir:ir+fieldsPerData-1, 1],
                                    DD[ir:ir+fieldsPerData-1, 2] )  # 4x6 matrix

        irow = (j-1)*nrow + 1
        icol = (j-1)*ncol + 1

        dZdD[irow:irow+nrow-1, icol:icol+ncol-1] = deriv1

        icol += nr*ncol
        dZdD[irow:irow+nrow-1, icol:icol+ncol-1] = deriv2

        j += 1
    end  # ir

    dropzeros!(dZdD)
    return dZdD
end # function calcFullZTEMDerivative

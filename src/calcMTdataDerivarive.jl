
export calcMTdata, calcFullMTDerivative

#--------------------------------------------------------------------------

function calcMTdata( DD::Array{Complex128,2} )  # (ndata,2) E and H data from polarizations 1 and 2
    # Data must be ordered [ ex,ey,hx,hy, ... ]
    # Calculate: Z = E H^-1

    #             ( Ex1 Ex2 )    ( Ex1*Hy2 - Ex2*Hy1   -Ex1*Hx2 + Ex2*Hx1 )
    # (Z11 Z12)   ( Ey1 Ey2 )    ( Ey1*Hy2 - Ey2*Hy1   -Ey1*Hx2 + Ey2*Hx1 )
    # (Z21 Z22) = -----------  =  -----------------------------------------
    #             ( Hx1 Hx2 )                Hx1*Hy2 - Hx2*Hy1
    #             ( Hy1 Hy2 )


    fieldsPerData = 4

    ndata = size(DD,1)   
    nr = div(ndata, fieldsPerData)

    zdpred = Array{Complex128}(2,2,nr)

    id = 0

    for ir = 1:fieldsPerData:ndata
        Ex1 = DD[ir,   1]
        Ey1 = DD[ir+1, 1]
        Hx1 = DD[ir+2, 1]
        Hy1 = DD[ir+3, 1]

        Ex2 = DD[ir,   2]
        Ey2 = DD[ir+1, 2]
        Hx2 = DD[ir+2, 2]
        Hy2 = DD[ir+3, 2]

        HH = Hx1*Hy2 - Hx2*Hy1

        id += 1
        zdpred[1,1,id] = (  Ex1*Hy2 - Ex2*Hy1 ) / HH
        zdpred[2,1,id] = (  Ey1*Hy2 - Ey2*Hy1 ) / HH
        zdpred[1,2,id] = ( -Ex1*Hx2 + Ex2*Hx1 ) / HH
        zdpred[2,2,id] = ( -Ey1*Hx2 + Ey2*Hx1 ) / HH
    end  # ir

    return zdpred
end  # function calcMTdata

#--------------------------------------------------------------------------

function calcFullMTDerivative( DD::Array{Complex128,2} )  # (nr*4,2) E and H data from polarizations 1 and 2
    # Data must be ordered [ ex,ey,hx,hy, ... ]
    # Calculate: Z = E H^-1

    # Output is a 2-block diagonal matrix.
    # In the derivative matrix, for each receiver, the 8 rows are ordered:
    #    [ Z11r Z11i  Z21r Z21i  Z12r Z12i  Z22r Z22i ]
    # and the 2*8 columns are ordered:
    #    [ Ex1r Ex1i  Ey1r Ey1i  Hx1r Hx1i  Hy1r Hy1i ]
    #    [ Ex2r Ex2i  Ey2r Ey2i  Hx2r Hx2i  Hy2r Hy2i ]

    fieldsPerData = 4

    ndata = size(DD,1) 
    if mod(ndata,fieldsPerData) != 0
       error("# of data should be divisible by 4.")
    end

    nr = div(ndata, fieldsPerData)  # number of data locations

    nrow =  8
    ncol =  8 #16
    dZdD = spzeros(Float64, nr*nrow, nr*ncol*2)

    j = 1
    for ir = 1:fieldsPerData:ndata

        deriv1,deriv2 = MTderivs( DD[ir:ir+fieldsPerData-1, 1],
                                 DD[ir:ir+fieldsPerData-1, 2] )  # 8x16 matrix

        irow = (j-1)*nrow + 1
        icol = (j-1)*ncol + 1

        dZdD[irow:irow+nrow-1, icol:icol+ncol-1] = deriv1

        icol += nr*ncol
        dZdD[irow:irow+nrow-1, icol:icol+ncol-1] = deriv2

        j += 1
    end  # ir

    dropzeros!(dZdD)
    return dZdD
end # function calcFullMTDerivative

#--------------------------------------------------------------------------

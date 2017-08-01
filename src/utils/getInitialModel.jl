export getInitialmodel

function getInitialmodel(M,itopo, halfSpaceCond, backCond = 1e-8)

	isactive = findActiveCells!(M,itopo)
	
	nactive    = sum(isactive)
	nnotactive = length(isactive) - nactive

	# model a half space
	sigma    = ones(nactive)*halfSpaceCond
	sigmaBck = ones(nnotactive)*backCond

	return sigma, sigmaBck, isactive

end
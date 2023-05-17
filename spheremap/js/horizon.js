
function eqToHorizon(ra, decl, lat, lst) {
    // Stupid simple rough approximation, ignores aberration, precession, diffraction, etc.
    const ha = lst - ra
    const alt = Math.asin(
        Math.sin(decl) * Math.sin(lat) + Math.cos(decl) * Math.cos(lat) * Math.cos(ha)
    )
    const az = Math.atan2(
        -1 * Math.cos(decl) * Math.cos(lat) * Math.sin(ha),
        Math.sin(decl) - Math.sin(lat) * Math.sin(alt)
    )
    const coords = [alt, az]
    return coords
}

function eqToHorizonCart(ra, decl, lat, lst) {
    const horizon = eqToHorizon(ra, decl, lat, lst)
    const alt = horizon[0]
    const az = horizon[1]
    const zd = Math.PI / 2 - alt
    const x = -1 * zd * Math.sin(az)
    const y = zd * Math.cos(az)
    let coords = [x, y]
    if ((x ** 2 + y ** 2) > ((Math.PI / 2) ** 2)) {
        coords[0] = NaN
        coords[1] = NaN
    }
    return coords
}

function horizonTransform() {
    const result = new Array(xs.length)

    if (xs.length == 0) {
        return result
    }

    const coord_idx = (proj_coord == 'x') ? 0 : 1

    
    // Columns can contain lists of coords (eg for corners of healpixels), or individual points.
    // If they are lists, iteratate over each element. Otherwise, just apply the rotation to the point.
    for (let i = 0; i < ra.length; i++) {
        if (typeof (ra[i]) === 'number') {
            const horizonCart = eqToHorizonCart(ra[i] * Math.PI / 180, decl[i] * Math.PI / 180, lat, lst)
            result[i] = horizonCart[coord_idx]
        } else {
            result[i] = new Array(ra[i].length)
            for (let j = 0; j < ra[i].length; j++) {
                const horizonCart = eqToHorizonCart(ra[i][j] * Math.PI / 180, decl[i][j] * Math.PI / 180, lat, lst)
                result[i][j] = horizonCart[coord_idx]
            }
        }
    }

    return result
}

function updateHorizonData() {
    console.log("Updating horizon data")
    for (let i = 0; i < ra.length; i++) {
        if (typeof (ra[i]) === 'number') {
            const horizonCart = eqToHorizonCart(ra[i] * Math.PI / 180, decl[i] * Math.PI / 180, lat, lst)
            if ('x_hz' in data) { data['x_hz'][i] = horizonCart[0] }
            if ('y_hz' in data) { data['y_hz'][i] = horizonCart[1] }
        } else {
            for (let j = 0; j < ra[i].length; j++) {
                const horizonCart = eqToHorizonCart(ra[i][j] * Math.PI / 180, decl[i][j] * Math.PI / 180, lat, lst)
                if ('x_hz' in data) { data['x_hz'][i][j] = horizonCart[0] }
                if ('y_hz' in data) { data['y_hz'][i][j] = horizonCart[1] }
            }
        }
    }

}

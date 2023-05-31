function eqToLambertAEA(ra, decl, hemisphere, west_right) {
    // Follow notation of Snyder p. 87-88
    let theta = (west_right === (hemisphere === 'south')) ? -1 * ra : ra
    const phi = (hemisphere === 'south') ? -1 * decl : decl

    // Choose an R to match that used by healpy
    const R = 1

    const rho = 2 * R * Math.sin((Math.PI / 2 - phi) / 2)
    let x = rho * Math.sin(theta)
    let y = rho * Math.cos(theta)
    if (phi > 89.9 * Math.PI / 180) {
        x = Math.NaN
        y = Math.NaN
    }
    return [x, y]
}

const hemisphere = (lat > 0) ? 'north' : 'south'

function laeaTransform() {
    const result = new Array(xs.length)

    if (xs.length == 0) {
        return result
    }
    
    const coord_idx = (proj_coord == 'x') ? 0 : 1

    for (let i = 0; i < ra.length; i++) {
        if (typeof (ra[i]) === 'number') {
            const laea = eqToLambertAEA(ra[i] * Math.PI / 180, decl[i] * Math.PI / 180, hemisphere, true)
            result[i] = laea[coord_idx]
        } else {
            result[i] = new Array(ra[i].length)
            for (let j = 0; j < ra[i].length; j++) {
                const laea = eqToLambertAEA(ra[i][j] * Math.PI / 180, decl[i][j] * Math.PI / 180, hemisphere, true)
                result[i][j] = laea[coord_idx]
            }
        }
    }

    return result
}

function updateLAEAData() {
    for (let i = 0; i < ra.length; i++) {
        if (typeof (ra[i]) === 'number') {
            const laea = eqToLambertAEA(ra[i] * Math.PI / 180, decl[i] * Math.PI / 180, hemisphere, true)
            if ('x_laea' in data) { data['x_laea'][i] = laea[0] }
            if ('y_laea' in data) { data['y_laea'][i] = laea[1] }
        } else {
            for (let j = 0; j < ra[i].length; j++) {
                const laea = eqToLambertAEA(ra[i][j] * Math.PI / 180, decl[i][j] * Math.PI / 180, hemisphere, true)
                if ('x_laea' in data) { data['x_laea'][i][j] = laea[0] }
                if ('y_laea' in data) { data['y_laea'][i][j] = laea[1] }
            }
        }
    }

}
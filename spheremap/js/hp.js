function eqToCart(ra, decl) {
    const theta = Math.PI / 2 - decl
    const phi = ra
    const x = Math.sin(theta) * Math.cos(phi)
    const y = Math.sin(theta) * Math.sin(phi)
    const z = Math.cos(theta)
    return [x, y, z]
}

function hpTransform() {

    if (xs.length == 0) {
        return new Array(xs.length)
    }
    
    const x_hp = new Array(ra.length)
    const y_hp = new Array(ra.length)
    const z_hp = new Array(ra.length)

    for (let i = 0; i < ra.length; i++) {
        if (typeof (ra[0]) === 'number') {
            let cartCoords = eqToCart(ra[i]*Math.PI/180, decl[i]*Math.PI/180)
            if (cartCoords[2]<0) {
                x_hp[i] = cartCoords[0]
                y_hp[i] = cartCoords[1]
                z_hp[i] = cartCoords[2]
            } else {
                x_hp[i] = NaN
                y_hp[i] = NaN
                z_hp[i] = NaN
            }
        } else {
            x_hp[i] = new Array(ra[i].length)
            y_hp[i] = new Array(ra[i].length)
            z_hp[i] = new Array(ra[i].length)
            for (let j = 0; j < ra[i].length; j++) {
                let cartCoords = eqToCart(ra[i][j]*Math.PI/180, decl[i][j]*Math.PI/180)
                if (cartCoords[2]<0) {
                    x_hp[i][j] = cartCoords[0]
                    y_hp[i][j] = cartCoords[1]
                    z_hp[i][j] = cartCoords[2]
                } else {
                    x_hp[i][j] = NaN
                    y_hp[i][j] = NaN
                    z_hp[i][j] = NaN
                }
            }
        }
    }

    const result = (proj_coord == 'x') ? x_hp : y_hp

    return result
}

function updateHpData() {

    let x_hp = NaN
    let y_hp = NaN
    let z_hp = NaN

    if ('x_hp' in data) {
        x_hp = data['x_hp']
    } else {
        x_hp = new Array(ra.length)
        for (let i = 0; i < ra.length; i++) {
            if (! typeof (ra[i]) === 'number') {
                x_hp[i] = new Array(ra[i].length)
            }
        }
    }

    if ('y_hp' in data) {
        y_hp = data['y_hp']
    } else {
        y_hp = new Array(ra.length)
        for (let i = 0; i < ra.length; i++) {
            if (! typeof (ra[i]) === 'number') {
                y_hp[i] = new Array(ra[i].length)
            }
        }
    }

    if ('z_hp' in data) {
        z_hp = data['z_hp']
    } else {
        z_hp = new Array(ra.length)
        for (let i = 0; i < ra.length; i++) {
            if (! typeof (ra[i]) === 'number') {
                z_hp[i] = new Array(ra[i].length)
            }
        }
    }

    for (let i = 0; i < ra.length; i++) {
        if (typeof (ra[i]) === 'number') {
            let cartCoords = eqToCart(ra[i]*Math.PI/180, decl[i]*Math.PI/180)
            if (cartCoords[2]<0) {
                x_hp[i] = cartCoords[0]
                y_hp[i] = cartCoords[1]
                z_hp[i] = cartCoords[2]
            } else {
                x_hp[i] = NaN
                y_hp[i] = NaN
                z_hp[i] = NaN
            }
        } else {
            for (let j = 0; j < ra[i].length; j++) {
                let cartCoords = eqToCart(ra[i][j]*Math.PI/180, decl[i][j]*Math.PI/180)
                if (cartCoords[2]<0) {
                    x_hp[i][j] = cartCoords[0]
                    y_hp[i][j] = cartCoords[1]
                    z_hp[i][j] = cartCoords[2]
                } else {
                    x_hp[i][j] = NaN
                    y_hp[i][j] = NaN
                    z_hp[i][j] = NaN
                }
            }
        }
    }

    // Columns can contain lists of coords (eg for corners of healpixels), or individual points.
    // If they are lists, iteratate over each element. Otherwise, just apply the rotation to the point.
    if (typeof (x_hp[0]) === 'number') {
        for (let i = 0; i < x_hp.length; i++) {
            if ('x_hp' in data) { data['x_hp'][i] = x_hp[i] }
            if ('y_hp' in data) { data['y_hp'][i] = y_hp[i] }
            if ('z_hp' in data) { data['z_hp'][i] = z_hp[i] }
        }
    } else {
        for (let i = 0; i < x_hp.length; i++) {
            for (let j = 0; j < x_hp[i].length; j++) {
                if ('x_hp' in data) { data['x_hp'][i][j] = x_hp[i][j] }
                if ('y_hp' in data) { data['y_hp'][i][j] = y_hp[i][j] }
                if ('z_hp' in data) { data['z_hp'][i][j] = z_hp[i][j] }
            }
        }
    }

}

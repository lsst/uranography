function eqToMollweide(ra, decl, west_right) {
    const tolerance = 0.001
    const max_iter = 1000
    const R = 1 / Math.sqrt(2)
    const dir_sign = west_right ? -1 : 1
    const wra = (ra + Math.PI) % (2 * Math.PI) - Math.PI

    // Return NaNs if near the discontinuity
    if (Math.abs(ra - Math.PI) < (Math.PI / 180) / Math.cos(decl)) {
        let xy = [Math.NaN, Math.NaN]
        return xy
    }

    function compute_xy(theta) {
        const x = dir_sign * R * (2 / Math.PI) * Math.sqrt(2) * wra * Math.cos(theta)
        const y = Math.sqrt(2) * R * Math.sin(theta)
        return [x, y]
    }

    let theta = decl
    let xy = compute_xy(theta)

    for (let iter = 1; iter < max_iter; iter++) {
        if (Math.cos(theta) ** 2 <= Math.cos(Math.PI / 2) ** 2) {
            // We are too close to a pole to refine further
            break
        }

        const theta0 = theta
        const xy0 = compute_xy(theta0)

        const delta_theta = -(2 * theta0 + Math.sin(2 * theta0) - Math.PI * Math.sin(decl)) / (
            4 * (Math.cos(theta0) ** 2)
        )
        theta = theta0 + delta_theta
        xy = compute_xy(theta)

        if ((Math.abs(xy[0] - xy0[0]) < tolerance) & (Math.abs(xy[1] - xy0[1]) < tolerance)) {
            break
        }
    }

    return xy
}

function mollweideTransform() {
    const result = new Array(xs.length)

    if (xs.length == 0) {
        return result
    }

    const coord_idx = (proj_coord == 'x') ? 0 : 1

    for (let i = 0; i < ra.length; i++) {
        if (typeof (ra[i]) === 'number') {
            const moll = eqToMollweide(ra[i] * Math.PI / 180, decl[i] * Math.PI / 180, true)
            result[i] = moll[coord_idx]
        } else {
            result[i] = new Array(ra[i].length)
            for (let j = 0; j < ra[i].length; j++) {
                const moll = eqToMollweide(ra[i][j] * Math.PI / 180, decl[i][j] * Math.PI / 180, true)
                result[i][j] = moll[coord_idx]
            }
        }
    }

    return result
}

function updateMollweideData() {

    for (let i = 0; i < ra.length; i++) {
        if (typeof (ra[i]) === 'number') {
            const moll = eqToMollweide(ra[i] * Math.PI / 180, decl[i] * Math.PI / 180, true)
            if ('x_moll' in data) { data['x_moll'][i] = moll[0] }
            if ('y_moll' in data) { data['y_moll'][i] = moll[1] }
        } else {
            for (let j = 0; j < ra[i].length; j++) {
                const moll = eqToMollweide(ra[i][j] * Math.PI / 180, decl[i][j] * Math.PI / 180, true)
                if ('x_moll' in data) { data['x_moll'][i][j] = moll[0] }
                if ('y_moll' in data) { data['y_moll'][i][j] = moll[1] }
            }
        }
    }

}
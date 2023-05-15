function get_gpu() {
    try {
        return gpu;
    } catch (e) {
        return new GPU();
    }
}

function rotateCart(ux, uy, uz, angle, x0, y0, z0) {
    // ux, uy, uz is a vector that defines the axis of rotation
    // angle is in radians
    // following https://en.wikipedia.org/wiki/Rotation_matrix

    const cosa = Math.cos(angle)
    const ccosa = 1 - cosa
    const sina = Math.sin(angle)
    const rxx = cosa + ux * ux * ccosa
    const rxy = ux * uy * ccosa - uz * sina
    const rxz = ux * uz * ccosa + uy * sina
    const ryx = uy * ux * ccosa + uz * sina
    const ryy = cosa + uy * uy * ccosa
    const ryz = uy * uz * ccosa - ux * sina
    const rzx = uz * ux * ccosa - uy * sina
    const rzy = uz * uy * ccosa + ux * sina
    const rzz = cosa + uz * uz * ccosa
    const x = rxx * x0 + rxy * y0 + rxz * z0
    const y = ryx * x0 + ryy * y0 + ryz * z0
    const z = rzx * x0 + rzy * y0 + rzz * z0
    return [x, y, z]
}

function applyRotations(hpx, hpy, hpz, codecl, ra, orient, npoleCoords1) {
    // We are looking out of the sphere from the inside, so the center is 180 degrees 
    // from the front of the sphere, hence the pi.
    const decl_rot = Math.PI + codecl
    const ra_rot = ra - Math.PI / 2

    const coords1 = rotateCart(1, 0, 0, decl_rot, hpx, hpy, hpz)
    const coords2 = rotateCart(npoleCoords1[0], npoleCoords1[1], npoleCoords1[2], ra_rot, coords1[0], coords1[1], coords1[2])
    let coords = rotateCart(0, 0, 1, orient, coords2[0], coords2[1], coords2[2])

    // In astronomy, we are looking out of the sphere from the center to the back
    // (which naturally results in west to the right).
    // Positive z is out of the screen behind us, and we are at the center,
    // so to visible part is when z is negative (coords[2]<=0).
    // So, stuff the points with positive z to NaN so they are
    // not shown, because they are behind the observer.

    // Use 5*Number.EPSILON instead of exactly 0, because the
    // assorted trig operations result in values slightly above or below
    // 0 when the horizon is in principle exactly 0, and this gives an
    // irregularly dotted/dashed appearance to the horizon if 
    // a cutoff of exactly 0 is used.
    if (coords[2] > 5 * Number.EPSILON) {
        coords[0] = NaN
        coords[1] = NaN
    }
    return coords
}

function horizonToEq(lat, alt, az, lst) {
    // Stupid simple rough approximation, ignores aberration, precession, diffraction, etc.
    // Doing this "correctly" would make this much more complicated and much slower, and
    // of the dates of relevance won't make a significant difference.
    const decl = Math.asin(Math.sin(alt) * Math.sin(lat) + Math.cos(lat) * Math.cos(alt) * Math.cos(az))
    const ha = Math.atan2(
        -1 * Math.cos(alt) * Math.cos(lat) * Math.sin(az),
        Math.sin(alt) - Math.sin(lat) * Math.sin(decl)
    )
    const ra = lst - ha
    const coords = [ra, decl]
    return coords
}

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

function eqToCart(ra, decl) {
    const theta = Math.PI / 2 - decl
    const phi = ra
    const x = Math.sin(theta) * Math.cos(phi)
    const y = Math.sin(theta) * Math.sin(phi)
    const z = Math.cos(theta)
    return [x, y, z]
}

function cartToEq(x, y, z) {
    const theta = Math.acos(z)
    const ra = Math.atan2(y, x)
    const decl = Math.PI / 2 - theta
    return [ra, decl]
}

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

function computeLocalSiderealTime(mjd, longitude) {
    // Computes the Mean Sidereal Time

    // Follow Meeus's _Astronomical_Algorithms_ 2nd edition, equation 12.4
    // Avoid obvious simplifications to make it easier to check the
    // numbers in the equation exactly.
    // Meeus follows tho IAU recommendation of 1982.
    const jd = mjd + 2400000.5
    const t = (jd - 2451545.0) / 36525.0
    const theta0 = 280.46061837 + (360.98564736629 * (jd - 2451545.0)) + (0.000387933 * t * t) - (t * t * t / 38710000)
    const lon_deg = longitude * (180.0 / Math.PI)
    const lst_deg = ((theta0 + lon_deg) % 360 + 360) % 360
    const lst = lst_deg * Math.PI / 180.0
    return lst
}

const multiplyMultiMatrix = get_gpu().createKernel(function (coeff, data) {
    let sum = 0;
    let hpix = this.thread.y
    let corner = this.thread.x
    let out_coord = this.thread.z
    for (let in_coord = 0; in_coord < 3; in_coord++) {
        sum += coeff[out_coord][in_coord] * data[in_coord][hpix][corner];
    }
    return sum;
})

function multiplyRotationMatrix(a, b) {
    let result = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    for (let x = 0; x < 3; x++) {
        for (let y = 0; y < 3; y++) {
            for (let i = 0; i < 3; i++) {
                result[y][x] = result[y][x] + a[i][x] * b[y][i]
            }
        }
    }
    return result;
}

function computeRotationMatrix(ux, uy, uz, angle) {
    const cosa = Math.cos(angle)
    const ccosa = 1 - cosa
    const sina = Math.sin(angle)
    const rxx = cosa + ux * ux * ccosa
    const rxy = ux * uy * ccosa - uz * sina
    const rxz = ux * uz * ccosa + uy * sina
    const ryx = uy * ux * ccosa + uz * sina
    const ryy = cosa + uy * uy * ccosa
    const ryz = uy * uz * ccosa - ux * sina
    const rzx = uz * ux * ccosa - uy * sina
    const rzy = uz * uy * ccosa + ux * sina
    const rzz = cosa + uz * uz * ccosa
    const matrix = [[rxx, rxy, rxz], [ryx, ryy, ryz], [rzx, rzy, rzz]]
    return matrix
}

function applyHealpixRotations(hpx, hpy, hpz, codecl, ra, orient, npoleCoords1) {
    // We are looking out of the sphere from the inside, so the center is 180 degrees 
    // from the front of the sphere, hence the pi.
    const decl_rot = Math.PI + codecl
    const ra_rot = ra - Math.PI / 2

    const matrix1 = computeRotationMatrix(1, 0, 0, decl_rot)
    const matrix2 = computeRotationMatrix(npoleCoords1[0], npoleCoords1[1], npoleCoords1[2], ra_rot)
    const matrix3 = computeRotationMatrix(0, 0, 1, orient)
    const rotMatrix = multiplyRotationMatrix(multiplyRotationMatrix(matrix1, matrix2), matrix3)
    let coords = multiplyMultiMatrix(rotMatrix, [hpx, hpy, hpz])

    // In astronomy, we are looking out of the sphere from the center to the back
    // (which naturally results in west to the right).
    // Positive z is out of the screen behind us, and we are at the center,
    // so to visible part is when z is negative (coords[2]<=0).
    // So, stuff the points with positive z to NaN so they are
    // not shown, because they are behind the observer.

    // Use 5*Number.EPSILON instead of exactly 0, because the
    // assorted trig operations result in values slightly above or below
    // 0 when the horizon is in principle exactly 0, and this gives an
    // irregularly dotted/dashed appearance to the horizon if 
    // a cutoff of exactly 0 is used.
    for (let hpix = 0; hpix < coords[0].length; hpix++) {
        for (let corner = 0; corner < coords[0][0].length; corner++) {
            if (coords[2][hpix][corner] > 5 * Number.EPSILON) {
                coords[0][hpix][corner] = NaN
                coords[1][hpix][corner] = NaN
            }
        }
    }
    return coords
}

const data = data_source.data

lat = lat * Math.PI / 180
lon = lon * Math.PI / 180

const alt = center_alt_slider.value * Math.PI / 180
const az = center_az_slider.value * Math.PI / 180
const mjd = mjd_slider.value
const lst = computeLocalSiderealTime(mjd, lon)

const eqCoords = horizonToEq(lat, alt, az, lst)
const ra = eqCoords[0]
const decl = eqCoords[1]
const codecl = Math.PI / 2 - decl
const npoleCoords1 = rotateCart(1, 0, 0, codecl, 0, 0, 1)

const hemisphere = (lat > 0) ? 'north' : 'south'

/* To get the orientation
   - Get the coordinates of a point we want to be directly "up" from center (slightly higher alt, same az)
   - Look where it would end up with orientation 0
   - Get the angle of the desired top with the actual top
   - Reverse to get the rotation
   */
let upAlt = alt + Math.PI / 180
let upAz = az
const upEq = horizonToEq(lat, upAlt, upAz, lst)
const upCart0 = eqToCart(upEq[0], upEq[1])
const upCart3 = applyRotations(upCart0[0], upCart0[1], upCart0[2], codecl, ra, 0, npoleCoords1)
const orient = Math.PI / 2 - Math.atan2(upCart3[1], upCart3[0])

function updateData() {
    let pointEq = NaN
    let eqUpdated = false
    let hzUpdated = true

    if (data['x_hp'].length == 0) {
        return
    }

    // Columns can contain lists of coords (eg for corners of healpixels), or individual points.
    // If they are lists, iteratate over each element. Otherwise, just apply the rotation to the point.
    if (typeof (data['x_hp'][0]) === 'number') {
        if ('alt' in data) {
            hzUpdated = false
            for (let i = 0; i < data['x_hp'].length; i++) {
                pointEq = horizonToEq(lat, data['alt'][i] * Math.PI / 180, data['az'][i] * Math.PI / 180, lst)
                let new_ra = pointEq[0] * 180 / Math.PI
                let new_decl = pointEq[1] * 180 / Math.PI
                // If the point does not actually change, be user eqUpdated remains false
                // so we do not call the expensive eqToMollweide below.
                if ((data['ra'][i] != new_ra) || (data['decl'][i] != new_decl)) {
                    data['ra'][i] = new_ra
                    data['decl'][i] = new_decl
                    eqUpdated = true

                    let cartCoords = eqToCart(pointEq[0], pointEq[1])
                    data['x_hp'][i] = cartCoords[0]
                    data['y_hp'][i] = cartCoords[1]
                    data['z_hp'][i] = cartCoords[2]
                }
            }
        }
        for (let i = 0; i < data['x_hp'].length; i++) {
            const coords = applyRotations(data['x_hp'][i], data['y_hp'][i], data['z_hp'][i], codecl, ra, orient, npoleCoords1)
            data['x_orth'][i] = coords[0]
            data['y_orth'][i] = coords[1]
            data['z_orth'][i] = coords[2]
        }

        if (eqUpdated && ('x_laea' in data)) {
            for (let i = 0; i < data['x_hp'].length; i++) {
                const laea = eqToLambertAEA(data['ra'][i] * Math.PI / 180, data['decl'][i] * Math.PI / 180, hemisphere, true)
                data['x_laea'][i] = laea[0]
                data['y_laea'][i] = laea[1]
            }
        }

        if (eqUpdated && ('x_moll' in data)) {
            for (let i = 0; i < data['x_hp'].length; i++) {
                const moll = eqToMollweide(data['ra'][i] * Math.PI / 180, data['decl'][i] * Math.PI / 180, true)
                data['x_moll'][i] = moll[0]
                data['y_moll'][i] = moll[1]
            }
        }

        if (hzUpdated && ('x_hz' in data)) {
            for (let i = 0; i < data['x_hp'].length; i++) {
                const horizonCart = eqToHorizonCart(data['ra'][i] * Math.PI / 180, data['decl'][i] * Math.PI / 180, lat, lst)
                data['x_hz'][i] = horizonCart[0]
                data['y_hz'][i] = horizonCart[1]
            }
        }
    } else {
        multiplyMultiMatrix.setOutput([data_source.data['x_hp'][0].length, data_source.data['x_hp'].length, 3]);

        if ('alt' in data) {
            hzUpdated = false
            for (let j = 0; j < data['x_hp'][0].length; j++) {
                for (let i = 0; i < data['x_hp'].length; i++) {
                    pointEq = horizonToEq(lat, data['alt'][i][j] * Math.PI / 180, data['az'][i][j] * Math.PI / 180, lst)
                    const new_ra = pointEq[0] * 180 / Math.PI
                    const new_decl = pointEq[1] * 180 / Math.PI
                    if ((data['ra'][i][j] != new_ra) || (data['decl'][i][j] != new_decl)) {
                        data['ra'][i][j] = pointEq[0] * 180 / Math.PI
                        data['decl'][i][j] = pointEq[1] * 180 / Math.PI
                        eqUpdated = True

                        let cartCoords = eqToCart(pointEq[0], pointEq[1])
                        data['x_hp'][i][j] = cartCoords[0]
                        data['y_hp'][i][j] = cartCoords[1]
                        data['z_hp'][i][j] = cartCoords[2]
                    }
                }
            }
        }

        const coords = applyHealpixRotations(data['x_hp'], data['y_hp'], data['z_hp'], codecl, ra, orient, npoleCoords1)
        for (let j = 0; j < data['x_hp'][0].length; j++) {
            for (let i = 0; i < data['x_hp'].length; i++) {
                // const coords = applyRotations(data['x_hp'][i][j], data['y_hp'][i][j], data['z_hp'][i][j], codecl, ra, orient, npoleCoords1)
                data['x_orth'][i][j] = coords[0][i][j]
                data['y_orth'][i][j] = coords[1][i][j]
                data['z_orth'][i][j] = coords[2][i][j]
            }
        }

        if (eqUpdated && ('x_laea' in data)) {
            for (let j = 0; j < data['x_hp'][0].length; j++) {
                for (let i = 0; i < data['x_hp'].length; i++) {
                    const laea = eqToLambertAEA(data['ra'][i][j] * Math.PI / 180, data['decl'][i][j] * Math.PI / 180, hemisphere, true)
                    data['x_laea'][i][j] = laea[0]
                    data['y_laea'][i][j] = laea[1]
                }
            }
        }

        if (eqUpdated && ('x_moll' in data)) {
            for (let j = 0; j < data['x_hp'][0].length; j++) {
                for (let i = 0; i < data['x_hp'].length; i++) {
                    const moll = eqToMollweide(data['ra'][i][j] * Math.PI / 180, data['decl'][i][j] * Math.PI / 180, true)
                    data['x_moll'][i][j] = moll[0]
                    data['y_moll'][i][j] = moll[1]
                }
            }
        }

        if (hzUpdated && ('x_hz' in data)) {
            for (let j = 0; j < data['x_hp'][0].length; j++) {
                for (let i = 0; i < data['x_hp'].length; i++) {
                    const horizonCart = eqToHorizonCart(data['ra'][i][j] * Math.PI / 180, data['decl'][i][j] * Math.PI / 180, lat, lst)
                    data['x_hz'][i][j] = horizonCart[0]
                    data['y_hz'][i][j] = horizonCart[1]
                }
            }
        }
    }

    if ('in_mjd_window' in data) {
        for (let i = 0; i < data['x_hp'].length; i++) {
            data['in_mjd_window'][i] = 0.5
            if ('min_mjd' in data) {
                if (mjd < data['min_mjd'][i]) {
                    data['in_mjd_window'][i] = 0.0
                }
            }
            if ('max_mjd' in data) {
                if (mjd > data['max_mjd'][i]) {
                    data['in_mjd_window'][i] = 0.0
                }
            }
        }
    }

    if (('recent_mjd' in data) && ('min_mjd' in data)) {
        for (let i = 0; i < data['x_hp'].length; i++) {
            let recent_mjd = 1.0 - (mjd - data['min_mjd'][i]) / data['fade_scale'][i]
            if ((recent_mjd < 0) || (recent_mjd > 1)) {
                recent_mjd = 0.0
            }
            data['recent_mjd'][i] = recent_mjd
        }
    }

}

updateData()

data_source.change.emit()

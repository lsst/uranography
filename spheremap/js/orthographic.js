
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

function get_gpu() {
    if (!window.hasOwnProperty('gpu')) {
        window.gpu = new GPU()
    }
    return window.gpu;
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


const alt = center_alt_slider.value * Math.PI / 180
const az = center_az_slider.value * Math.PI / 180

const centerEqCoords = horizonToEq(lat, alt, az, lst)
const centerRA = centerEqCoords[0]
const centerDecl = centerEqCoords[1]
const centerCodecl = Math.PI / 2 - centerDecl
const npoleCoords1 = rotateCart(1, 0, 0, centerCodecl, 0, 0, 1)

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
const upCart3 = applyRotations(upCart0[0], upCart0[1], upCart0[2], centerCodecl, centerRA, 0, npoleCoords1)
const orient = Math.PI / 2 - Math.atan2(upCart3[1], upCart3[0])

function orthoTransform() {
    const result = new Array(xs.length)

    if (xs.length == 0) {
        return result
    }

    const coord_idx = (proj_coord == 'x') ? 0 : 1

    const x_hp = new Array(ra.length)
    const y_hp = new Array(ra.length)
    const z_hp = new Array(ra.length)

    for (let i = 0; i < ra.length; i++) {
        if (typeof (ra[0]) === 'number') {
            let cartCoords = eqToCart(ra[i] * Math.PI / 180, decl[i] * Math.PI / 180)
            x_hp[i] = cartCoords[0]
            y_hp[i] = cartCoords[1]
            z_hp[i] = cartCoords[2]
        } else {
            x_hp[i] = new Array(ra[i].length)
            y_hp[i] = new Array(ra[i].length)
            z_hp[i] = new Array(ra[i].length)
            for (let j = 0; j < ra[i].length; j++) {
                let cartCoords = eqToCart(ra[i][j] * Math.PI / 180, decl[i][j] * Math.PI / 180)
                x_hp[i][j] = cartCoords[0]
                y_hp[i][j] = cartCoords[1]
                z_hp[i][j] = cartCoords[2]
            }
        }
    }

    // Columns can contain lists of coords (eg for corners of healpixels), or individual points.
    // If they are lists, iteratate over each element. Otherwise, just apply the rotation to the point.
    if (typeof (x_hp[0]) === 'number') {
        for (let i = 0; i < x_hp.length; i++) {
            const coords = applyRotations(x_hp[i], y_hp[i], z_hp[i], centerCodecl, centerRA, orient, npoleCoords1)
            result[i] = coords[coord_idx]
        }
    } else {
        multiplyMultiMatrix.setOutput([x_hp[0].length, x_hp.length, 3]);
        const coords = applyHealpixRotations(x_hp, y_hp, z_hp, centerCodecl, centerRA, orient, npoleCoords1)
        for (let i = 0; i < x_hp.length; i++) {
            result[i] = new Array(x_hp[i].length)
            for (let j = 0; j < x_hp[i].length; j++) {
                result[i][j] = coords[coord_idx][i][j]
            }
        }
    }

    return result
}

function updateOrthoData() {

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
            let cartCoords = eqToCart(ra[i] * Math.PI / 180, decl[i] * Math.PI / 180)
            x_hp[i] = cartCoords[0]
            y_hp[i] = cartCoords[1]
            z_hp[i] = cartCoords[2]
        } else {
            for (let j = 0; j < ra[i].length; j++) {
                let cartCoords = eqToCart(ra[i][j] * Math.PI / 180, decl[i][j] * Math.PI / 180)
                x_hp[i][j] = cartCoords[0]
                y_hp[i][j] = cartCoords[1]
                z_hp[i][j] = cartCoords[2]
            }
        }
    }

    // Columns can contain lists of coords (eg for corners of healpixels), or individual points.
    // If they are lists, iteratate over each element. Otherwise, just apply the rotation to the point.
    if (typeof (x_hp[0]) === 'number') {
        for (let i = 0; i < x_hp.length; i++) {
            const coords = applyRotations(x_hp[i], y_hp[i], z_hp[i], centerCodecl, centerRA, orient, npoleCoords1)
            if ('x_orth' in data) { data['x_orth'][i] = coords[0] }
            if ('y_orth' in data) { data['y_orth'][i] = coords[1] }
            if ('z_orth' in data) { data['z_orth'][i] = coords[2] }
        }
    } else {
        multiplyMultiMatrix.setOutput([x_hp[0].length, x_hp.length, 3]);
        const coords = applyHealpixRotations(x_hp, y_hp, z_hp, centerCodecl, centerRA, orient, npoleCoords1)
        for (let i = 0; i < x_hp.length; i++) {
            for (let j = 0; j < x_hp[i].length; j++) {
                if ('x_orth' in data) { data['x_orth'][i][j] = coords[0][i][j] }
                if ('y_orth' in data) { data['y_orth'][i][j] = coords[1][i][j] }
                if ('z_orth' in data) { data['z_orth'][i][j] = coords[2][i][j] }
            }
        }
    }

}

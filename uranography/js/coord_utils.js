
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

let data = {}
try {
    data = data_source.data
} catch(e) {
    console.log("****************HERE!*****************")
    // If we were not passed a data_source with data, assume we were 
    // passed a column of coordinate pairs to transform, and build our
    // data from that.
    data['ra'] = new Array(xs.length)
    data['decl'] = new Array(xs.length)
    for(let i=0; i<xs.length; i++) {
        data['ra'][i] = xs[i][0]
        data['decl'][i] = xs[i][1]
    }
}

const mjd = mjd_slider.value

lat = lat * Math.PI / 180
lon = lon * Math.PI / 180

const lst = computeLocalSiderealTime(mjd, lon)
let pointEq = NaN
let ra = NaN
let decl = NaN

if ('alt' in data) {
    ra = new Array(data['alt'].length)
    decl = new Array(data['alt'].length)

    for (let i = 0; i < data['alt'].length; i++) {
        if (typeof (data['alt'][0]) === 'number') {
            pointEq = horizonToEq(lat, data['alt'][i] * Math.PI / 180, data['az'][i] * Math.PI / 180, lst)
            ra[i] = pointEq[0] * 180 / Math.PI
            decl[i] = pointEq[1] * 180 / Math.PI
        } else {
            ra[i] = new Array(data['alt'][i].length)
            decl[i] = new Array(data['alt'][i].length)
            for (let j = 0; j < data['alt'][i].length; j++) {
                pointEq = horizonToEq(lat, data['alt'][i] * Math.PI / 180, data['az'][i] * Math.PI / 180, lst)
                ra[i][j] = pointEq[0] * 180 / Math.PI
                decl[i][j] = pointEq[1] * 180 / Math.PI
            }
        }
    }
} else {
    ra = data['ra']
    decl = data['decl']
}

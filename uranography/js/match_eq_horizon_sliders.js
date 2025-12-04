// Bokeh callback content to coordinate equatorial and horizon coordinate sliders

// If we are in the middle of executing this callback, do not
// run the callback again for the sliders adjusted by the
// initial callback: the initial callback will do all the
// work.
if (update_guard.text === "true") {
    return;
}
update_guard.text = "true";


const mjd = mjd_slider.value
const lat = lat_deg * Math.PI/180
const lon = lon_deg * Math.PI/180

// Follow Meeus's formula for converting
// MJD and longitude to LST.
// Not super high precision, but good
// enough for this purpose.
const jd = mjd + 2400000.5
const t = (jd - 2451545.0) / 36525.0
const theta0 = (
    280.46061837
    + (360.98564736629 * (jd - 2451545.0))
    + (0.000387933 * t * t)
    - (t * t * t / 38710000)
    )
const lst_deg = ((theta0 + lon_deg) % 360 + 360) % 360
const lst = lst_deg * Math.PI / 180.0

if (coord_to_update === 'eq') {
    // User adjusted a horizon slider, so
    // this callback should adjust the
    // equatorial sliders to match
    const alt = alt_slider.value * Math.PI / 180
    const az = az_slider.value * Math.PI / 180
    const decl = Math.asin(
        Math.sin(alt) * Math.sin(lat)
        + Math.cos(lat) * Math.cos(alt) * Math.cos(az)
        )
    const ha = Math.atan2(
        -1 * Math.cos(alt) * Math.cos(lat) * Math.sin(az),
        Math.sin(alt) - Math.sin(lat) * Math.sin(decl)
    )
    const ra_not_norm = lst - ha
    let ra = Math.atan2(Math.sin(ra_not_norm), Math.cos(ra_not_norm))
    if (ra <0 ) {
        ra = ra + 2 * Math.PI
    }
    const ra_deg = ra * 180 / Math.PI
    const decl_deg = decl * 180 / Math.PI

    // Only update if change is significant to avoid recursive updates
    const dra = ra_deg - ra_slider.value
    const ddecl = decl_deg - decl_slider.value
    if (dra**2 >= 1) {
        ra_slider.value = ra_deg
    }
    if (ddecl**2 >= 1) {
        decl_slider.value =  decl_deg
    }
} else if (coord_to_update === 'horizon') {
    // User adjusted an equatorial slider, so
    // this callback should adjust the
    // horizon sliders to match
    const ra = ra_slider.value * Math.PI / 180
    const decl = decl_slider.value * Math.PI / 180
    const ha = lst - ra
    const alt = Math.asin(
        Math.sin(decl) * Math.sin(lat) + Math.cos(decl) * Math.cos(lat) * Math.cos(ha)
    )
    let az = Math.atan2(
        -1 * Math.cos(decl) * Math.cos(lat) * Math.sin(ha),
        Math.sin(decl) - Math.sin(lat) * Math.sin(alt)
    )
    if (az < 0) {
        az = az + 2 * Math.PI
    }
    const alt_deg = alt * 180 / Math.PI
    const az_deg = az * 180 / Math.PI

    // only update if the change is significant to avoid recursive updates
    const dalt = alt_deg - alt_slider.value
    const daz = az_deg - alt_slider.value
    if (dalt**2 > 1) {
        alt_slider.value = alt_deg
    }
    if (daz**2 > 1) {
        az_slider.value = az_deg
    }
} else {
    console.log("Invalide coord_to_update")
}

function clear_guard() {
    update_guard.text = "false"
}
// If we run clear_guard directly, gets called before bokeh's
// queued triggers get run, so the guard is ineffective. Instead,
// use setTimeout to put it in the browsers event loop
// where it will get called after the bokeh callbacks
// get executed.
setTimeout(clear_guard, 0);

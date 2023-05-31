function updateData() {
    if ('ra' in data) { data['ra'] = ra }
    if ('decl' in data) { data['decl'] = decl }

    const update_hp = ('x_hp' in data) || ('y_hp' in data) || ('z_hp' in data)
    const update_orth = ('x_orth' in data) || ('y_orth' in data) || ('z_orth' in data)
    if (update_hp || update_orth) {
        updateOrthoData()
    }

    if (('x_hz' in data) || ('y_hz' in data)) {
        updateHorizonData()
    }

    if (('x_moll' in data) || ('y_moll' in data)) {
        updateMollweideData()
    }

    if (('x_laea' in data) || ('y_laea' in data)) {
        updateLAEAData()
    }
}

updateData()

data_source.change.emit()

; Script to convert AO-formatted catalogs to FHD formats
; Written by Ruby Byrne, 12/18

pro catalog_reformat_ao_to_fhd

    filepath = '/Users/ruby/EoR/extended_source_models_from_Ben_Fall2018/CenA_01_spec_index.txt'

    flux_struct={flux,xx:0.,yy:0.,xy:Complex(0.),yx:Complex(0.),I:0.,Q:0.,U:0.,V:0.}
    shape_struct={shape,x:0.,y:0.,angle:0.} ;gets used if keyword gaussian_source_models is set
    ;flux order is 0-3: xx, yy, xy, yx in apparent brightness; 4-7: I, Q, U, V in sky brightness
    ;flag type codes are 0: no flag, 1: low confidence 2: sidelobe contamination
    struct_base={id:-1L,x:0.,y:0.,ra:0.,dec:0.,ston:0.,freq:180.,alpha:0.,gain:1.,flag:0,extend:Ptr_new(),flux:flux_struct,shape:shape_struct}

    nlines = file_lines(filepath)
    file_contents = strarr(nlines)
    openr, datafile, filepath, /get_lun
    readf, datafile, file_contents
    free_lun, datafile

    source_start_lines = []
    component_start_lines = []
    for file_line = 0, nlines-1 do begin
      if strpos(file_contents[file_line], 'source {') ne -1 then begin
        source_start_lines = [source_start_lines, file_line]
      endif
      if strpos(file_contents[file_line], 'component {') ne -1 then begin
        component_start_lines = [component_start_lines, file_line]
      endif
    endfor

    target_freq_MHz = 180. ;If more than one measurement exists, use the one taken at the freq closest to this value
    source_struct=Replicate(struct_base,n_elements(source_start_lines)>1)
    ; Iterate over sources
    for source=0,n_elements(source_start_lines)-1 do begin
        source_start = source_start_lines[source]
        if source eq n_elements(source_start_lines)-1 then source_end=nlines-1 else source_end=source_start_lines[source+1]

        source_comp_start_lines = []
        for comp=0,n_elements(component_start_lines)-1 do begin
            if component_start_lines[comp] gt source_start and component_start_lines[comp] lt source_end then source_comp_start_lines = [source_comp_start_lines,component_start_lines[comp]]
        endfor
        if n_elements(source_comp_start_lines) eq 0 then use_source_comp_start_lines = [source_start] else use_source_comp_start_lines = source_comp_start_lines

        comp_struct=Replicate(struct_base,n_elements(use_source_comp_start_lines)>1)
        ; Iterate over components
        for comp=0,n_elements(use_source_comp_start_lines)-1 do begin
            comp_start = use_source_comp_start_lines[comp]
            if comp eq n_elements(use_source_comp_start_lines)-1 then comp_end=source_end else comp_end=use_source_comp_start_lines[comp+1]

            shape = []
            position = []
            frequency = []
            fluxdensity = []
            for file_line=comp_start,comp_end do begin
                if strpos(file_contents[file_line], 'shape') ne -1 then shape = [shape, file_contents[file_line]]
                if strpos(file_contents[file_line], 'position') ne -1 then position = [position, file_contents[file_line]]
                if strpos(file_contents[file_line], 'frequency') ne -1 then frequency = [frequency, file_contents[file_line]]
                if strpos(file_contents[file_line], 'fluxdensity') ne -1 then fluxdensity = [fluxdensity, file_contents[file_line]]
            endfor

            if n_elements(shape) eq 1 then begin
                shape_params = strsplit(shape, /extract)
                comp_struct[comp].shape.x=double(shape_params[1])
                comp_struct[comp].shape.y=double(shape_params[2])
                comp_struct[comp].shape.angle=double(shape_params[3])
            endif else begin
                if n_elements(shape) gt 1 then begin
                    print, 'ERROR: Ambiguous shape parameters provided.'
                    return
                endif
            endelse

            if n_elements(position) eq 1 then begin
                pos_params = strsplit(position, /extract)
                ra_hours = fix((strsplit(pos_params[1], 'h', /extract))[0])
                ra_minutes = fix((strsplit((strsplit(pos_params[1], 'h', /extract))[1], 'm', /extract))[0])
                ra_seconds = float((strsplit((strsplit(pos_params[1], 'm', /extract))[1], 's', /extract))[0])
                if ra_hours lt 0 then begin
                  ra_minutes = -ra_minutes
                  ra_seconds = -ra_seconds
                endif
                ra = ra_hours*15. + ra_minutes*(.25) + ra_seconds*(.05/12.)
                dec_degrees = fix((strsplit(pos_params[2], 'd', /extract))[0])
                dec_minutes = fix((strsplit((strsplit(pos_params[2], 'd', /extract))[1], 'm', /extract))[0])
                dec_seconds = float((strsplit((strsplit(pos_params[2], 'm', /extract))[1], 's', /extract))[0])
                if dec_degrees lt 0 then begin
                  dec_minutes = -dec_minutes
                  dec_seconds = -dec_seconds
                endif
                dec = dec_degrees + dec_minutes/60. + dec_seconds/3600.
                comp_struct[comp].ra = ra
                comp_struct[comp].dec = dec
            endif else begin
                if n_elements(position) eq 0 then begin
                    print, 'ERROR: No position parameter provided.'
                    return
                endif $
                else if n_elements(position) gt 1 then begin
                    print, 'ERROR: Ambiguous position parameters provided.'
                    return
                endif
            endelse

            if n_elements(frequency) eq 1 then begin
                freq_params = strsplit(frequency, /extract)
                comp_struct[comp].freq=float(freq_params[1])
            endif else begin
                if n_elements(frequency) eq 0 then print, 'WARNING: No frequency parameter provided.' $
                else if n_elements(frequency) gt 1 then begin
                    print, 'ERROR: Ambiguous frequency parameters provided.'
                    return
                endif
            endelse

            if n_elements(fluxdensity) eq 1 then begin
                flux_params = strsplit(fluxdensity, /extract)
                comp_struct[comp].flux.I=double(flux_params[2])
            endif else begin
                if n_elements(fluxdensity) eq 0 then begin
                    print, 'ERROR: No flux density parameter provided.'
                    return
                endif $
                else if n_elements(fluxdensity) gt 1 then begin
                    print, 'ERROR: Ambiguous flux density parameters provided.'
                    return
                endif
            endelse

        ; End iterate over components
        endfor

        if n_elements(comp_struct) eq 1 then source_struct[source]=comp_struct[0] else begin
            source_struct[source].extend = ptr_new(comp_struct)
            flux_total=0.
            weighted_ra = 0.
            weighted_dec = 0.
            for comp=0,n_elements(comp_struct)-1 do begin
                flux_total+=comp_struct[comp].flux.I
                weighted_ra+=(comp_struct[comp].flux.I*comp_struct[comp].ra)
                weighted_dec+=(comp_struct[comp].flux.I*comp_struct[comp].dec)
            endfor
            source_struct[source].ra = weighted_ra/flux_total
            source_struct[source].dec = weighted_dec/flux_total
            source_struct[source].flux.I = flux_total
        endelse
        
    ; End iterate over sources
    endfor
    
    catalog = source_struct
    save, catalog, filename='/Users/ruby/EoR/extended_source_models_from_Ben_Fall2018/CenA_gaussian_model.sav'

end

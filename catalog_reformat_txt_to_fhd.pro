; Script to convert from pyuvsim .txt catalog files to FHD readable format
; Written by Dara Storer, April 2019
; 
; Expecting .txt file of format:
; SOURCE_ID     RA[deg]     DEC[deg]    Flux[Jy]    Frequency[Hz]
; Assumes header occupies the first line only

pro catalog_reformat_txt_to_fhd

    ;Path to catalog to be reformatted
    filepath = '/Users/dstorer/Files/FHD_Pyuvsim_ComparisonSimulation/mock_catalog_heratext_2458098.38824015.txt'

    flux_struct={flux,xx:0.,yy:0.,xy:Complex(0.),yx:Complex(0.),I:0.,Q:0.,U:0.,V:0.}
    struct_base = {id:'',ra:0.,dec:0.,freq:0.,alpha:0.,extend:Ptr_new(),flux:flux_struct} ;structure of output catalog
    
    nlines = file_lines(filepath)
    file_contents = strarr(nlines)
    openr, datafile, filepath, /get_lun
    readf, datafile, file_contents
    free_lun, datafile
    
    source_struct = Replicate(struct_base, (nlines)-1)

    ; Populate sources 
    for source=0, nlines-2 do begin
        line = file_contents[source+1]
        split = strsplit(line, /extract)
        source_struct[source].id = split[0]
        source_struct[source].ra = float(split[1])
        source_struct[source].dec = float(split[2])
        source_struct[source].flux.I = float(split[3])
        source_struct[source].freq = float(100)
        source_struct[source].extend = ptr_new() 
    endfor
    
    catalog = source_struct
    save, catalog, filename='/Users/dstorer/Files/FHD_Pyuvsim_ComparisonSimulation/mock_catalog_heratext_2458098_38824015_test.sav'

end
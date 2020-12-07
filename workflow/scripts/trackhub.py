from pathlib import PurePath
import colorsys

import trackhub

def get_color(i, n_colors):
    h =  i/float(n_colors)
    s = 0.5 # color_intensity[filetype]
    color = [255*v for v in colorsys.hsv_to_rgb(h, 0.5, 0.5)]
    return ",".join([str(int(c)) for c in list(color)])

def histone_track(name, color="200,0,0"):
    composite = trackhub.CompositeTrack(
        name=name+"_composite",
        short_label=name,
        tracktype='bigWig',
        visibility='full',
        color=color
    )
    signal_view = trackhub.ViewTrack(
        name=name+"_signal",
        view='signal',
        visibility='full',
        tracktype='bigWig',
        short_label=name+'_signal',
        autoScale="on"
    )
    regions_view = trackhub.ViewTrack(
        name=name+'_region',
        view='regions',
        visibility='dense',
        tracktype='bigWig',
        short_label=name+'_regions')

    composite.add_view(signal_view)
    composite.add_view(regions_view)

    for signal_type in ["treat_pileup", "control_lambda"]:
        track = trackhub.Track(
            tracktype='bigWig',
            name=name+"_"+signal_type,
            url="%s_%s.bw" % (name, signal_type),
            short_label=signal_type,
            autoScale="on"
        )
        signal_view.add_tracks(track)

    for region_type in ["domains"]:
        track = trackhub.Track(
            name=name+"_"+region_type,
            url="%s_%s.bb" %(name, region_type),
            short_label=region_type,
            tracktype='bigBed')
        regions_view.add_tracks(track)
    return composite

names = [PurePath(name).stem.replace("_domains", "") for name in snakemake.input.domains]
colors = [get_color(i, len(names)) for i in range(len(names))]
hub, genomes_file, genome, trackdb = trackhub.default_hub(
    hub_name="testing",
    genome="hg38",
    email="knutdrand@gmail.com")
trackdb.add_tracks([histone_track(*pair) for pair in zip(names, colors)])
open(snakemake.output[0], "w").write(str(trackdb))

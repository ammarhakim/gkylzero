
# name="outputs/gk_wham_1x2v_p1_adiabatic"
# name="outputs/gk_mirror_adiabatic_elc_1x2v_p1_nonuniform_nosource"
name="gk_wham_1x2v_p1_boltz"
species="ion"

# Animations of distribution functions with vpar on z and mu on z
# pgkyl "$name-ion_[0-9]*.gkyl"\
#   interp -b gkhyb -p1 integrate 2 ev 'f[:] abs' animate --logz --zmin 1e-20 --fps 4 \
#   --saveas "$name vpar.mp4" &
# pgkyl "$name-ion_[0-9]*.gkyl"\
#   interp -b gkhyb -p1 integrate 1 ev 'f[:] abs' animate --logz --zmin 1e-4 --fps 4\
#   --saveas "$name mu.mp4" &

# Plot single frames of the distribution function
# pgkyl "$name-ion_71.gkyl" interp -b gkhyb -p1 integrate 2 ev 'f[:] abs' sel --z0 0.0 pl --logy &
frame=0
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 integrate 2 ev 'f[:] abs' pl --logz --zmin 1e-20 --title "frame $frame vpar"&
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 integrate 1 ev 'f[:] abs' pl --logz --zmin 1e-4 --title "frame $frame mu"&
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 sel --z0 0.0 ev 'f[:] abs' pl --logz --zmin 1e-10 --title "frame $frame z=0"&
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 sel --z0 0.90 ev 'f[:] abs' pl --logz --zmin 1e-10 --title "frame $frame z=.9"&
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 sel --z0 1.0 ev 'f[:] abs' pl --logz --zmin 1e-10 --title "frame $frame z=1.0"&
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 integrate 2 pl --logz --zmin 1e-20 --title "frame $frame vpar"&
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 integrate 1 pl --logz --zmin 1e-4 --title "frame $frame mu"&
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 sel --z0 0.0 pl --logz --zmin 1e-10 --title "frame $frame z=0"&

# Plot geometry quantities
pgkyl "$name-jacobgeo.gkyl" interp -b ms -p1 pl --title "jacobgeo"&
pgkyl "$name-jacobtot.gkyl" interp -b ms -p1 pl --title "jacobtot"&
pgkyl "$name-jacobtot_inv.gkyl" interp -b ms -p1 pl --title "jacobtot_inv"&
pgkyl "$name-jacogeo_inv.gkyl" interp -b ms -p1 pl --title "jacobgeo_inv"&
pgkyl "$name-b_i.gkyl" interp -b ms -p1 pl --title "b_i"&
pgkyl "$name-mapc2p.gkyl" interp -b ms -p1 pl --title "mapc2p"&
pgkyl "$name-bmag.gkyl" interp -b ms -p1 pl --title "bmag"&
pgkyl "$name-bmag_inv.gkyl" interp -b ms -p1 pl --title "bmag_inv"&
pgkyl "$name-bmag_inv_sq.gkyl" interp -b ms -p1 pl --title "bmag_inv_sq"&
pgkyl "$name-cmag.gkyl" interp -b ms -p1 pl --title "cmag"&

# Plot distribution function integrals
# $name-"$species"_0.gkyl interp -b gkhyb -p1 integrate 2 ev 'f[:] abs'  pl --logz --zmin 1e-20 &
# $name-"$species"_4.gkyl interp -b gkhyb -p1 integrate 2 ev 'f[:] abs'  pl --logz --zmin 1e-20 --title '4 micros vpar' &
# $name-"$species"_4.gkyl interp -b gkhyb -p1 integrate 1 ev 'f[:] abs'  pl --logz --zmin 1e-4  --title '4 micros mu' &
# $name-"$species"_1.gkyl interp -b gkhyb -p1 integrate 2 ev 'f[:] abs'  pl --logz --zmin 1e-20 --title 'without souce 1 micros' &
# $name-"$species"_1.gkyl interp -b gkhyb -p1 integrate 1 ev 'f[:] abs'  pl --logz --zmin 1e-4  --title 'without source 1 micros' &
# $name-"$species"_19.gkyl interp -b gkhyb -p1 integrate 2 ev 'f[:] abs'  pl --logz --zmin 1e-20 --title 'without souce 19 micros' &
# $name-"$species"_19.gkyl interp -b gkhyb -p1 integrate 1 ev 'f[:] abs'  pl --logz --zmin 1e-4  --title 'without source 19 micros' &

# Plot prim-moms
# $name-"$species"_nu_prim_moms_0.gkyl interp -b ms -p1 pl --title '0' &
# "$name-"$species"_prim_moms_[0-9]*.gkyl" interp -b ms -p1 collect pl --title 'prim--moms' &

# Distribution function at a single velocity space point at all z
# "$name-"$species"_[0-9]*.gkyl"\
#   interp -b gkhyb -p1 sel --z1 5e6 --z2 3.5e-15 ev 'f[:] abs' collect pl --logz --zmin 1e-4 \

#Plot phi
frame=117
# pgkyl "$name-field_[0-9]*.gkyl" interp -b ms -p1 collect sel --z1 0 pl --title 'phi' &
# pgkyl "$name-field_[0-9]*.gkyl" interp -b ms -p1 collect pl --title 'phi' &
# pgkyl "$name-field_$frame.gkyl" interp -b ms -p1 ev 'f[:] 940 /' pl --title "phi at frame $frame in unite e phi/Te" &

# Moments of the distribution function
frame=0
if [ "$species" = "elc" ]; then mass=9.11e-31
elif [ "$species" = "ion" ]; then mass=3.34e-27
else echo "species must be ion or elc"
fi
# pgkyl $name-"$species"_M0_$frame.gkyl interp -b ms -p1 pl --title 'Density' &
# pgkyl $name-"$species"_M0_$frame.gkyl $name-"$species"_M1_$frame.gkyl interp -b ms -p1 ev "f[1] f[0] /" pl --title 'Upar' &
# pgkyl $name-"$species"_M0_$frame.gkyl $name-"$species"_M2perp_$frame.gkyl interp -b ms -p1 ev "$mass f[1] f[0] / * 1.6e-19 /" pl --title 'Tperp (eV)' &
# pgkyl $name-"$species"_M0_$frame.gkyl $name-"$species"_M1_$frame.gkyl $name-"$species"_M2par_$frame.gkyl interp -b ms -p1 ev "f[2] f[1] f[1] * f[0] / - $mass * f[0] / 1.6e-19 /" pl --title 'Tpar (eV)' &

# Moments of the distribution function
frame1=0
frame2=25
if [ "$species" = "elc" ]; then mass=9.11e-31
elif [ "$species" = "ion" ]; then mass=3.34e-27
else echo "species must be ion or elc"
fi
# # # Compute frame 1
# pgkyl $name-"$species"_M0_$frame1.gkyl $name-"$species"_M1_$frame1.gkyl interp -b ms -p1 ev "f[1] f[0] /"  writ -f "$name-"$species"_upar_$frame1.gkyl" &
# pgkyl $name-"$species"_M0_$frame1.gkyl $name-"$species"_M2perp_$frame1.gkyl interp -b ms -p1 ev "$mass f[1] f[0] / * 1.6e-19 /" writ -f "$name-"$species"_tperp_$frame1.gkyl" &
# pgkyl $name-"$species"_M0_$frame1.gkyl $name-"$species"_M1_$frame1.gkyl $name-"$species"_M2par_$frame1.gkyl interp -b ms -p1 ev "f[2] f[1] f[1] * f[0] / - $mass * f[0] / 1.6e-19 /" writ -f "$name-"$species"_tpar_$frame1.gkyl" &
# # Compute frame 2
# pgkyl $name-"$species"_M0_$frame2.gkyl $name-"$species"_M1_$frame2.gkyl interp -b ms -p1 ev "f[1] f[0] /"  writ -f "$name-"$species"_upar_$frame2.gkyl" &
# pgkyl $name-"$species"_M0_$frame2.gkyl $name-"$species"_M2perp_$frame2.gkyl interp -b ms -p1 ev "$mass f[1] f[0] / * 1.6e-19 /" writ -f "$name-"$species"_tperp_$frame2.gkyl" &
# pgkyl $name-"$species"_M0_$frame2.gkyl $name-"$species"_M1_$frame2.gkyl $name-"$species"_M2par_$frame2.gkyl interp -b ms -p1 ev "f[2] f[1] f[1] * f[0] / - $mass * f[0] / 1.6e-19 /" writ -f "$name-"$species"_tpar_$frame2.gkyl" &
# # Plot the results
# pgkyl $name-"$species"_upar_$frame1.gkyl $name-"$species"_upar_$frame2.gkyl pl --title 'Upar' -f0 &
# pgkyl $name-"$species"_tperp_$frame1.gkyl $name-"$species"_tperp_$frame2.gkyl pl --title 'Tperp (eV)' -f0 &
# pgkyl $name-"$species"_tpar_$frame1.gkyl $name-"$species"_tpar_$frame2.gkyl pl --title 'Tpar (eV)' -f0 &
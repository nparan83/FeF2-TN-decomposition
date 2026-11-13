# fef2_setup.sh - pseudopotentials for FeF2 (PSlibrary PAW). Run once on login node. ASCII only.
cd /arf/scratch/nparan/FeF2 && mkdir -p pseudo && cd pseudo
for f in Fe.pbe-spn-kjpaw_psl.0.2.1.UPF F.pbe-n-kjpaw_psl.0.1.UPF; do
  [ -s $f ] || cp /arf/scratch/nparan/lnmo/pseudo/$f . 2>/dev/null
  [ -s $f ] || curl -sfL -o $f https://pseudopotentials.quantum-espresso.org/upf_files/$f
  grep -q "<PP_HEADER" $f 2>/dev/null || { echo "indirilemedi: $f -> ld1.x ile uretilecek"; rm -f $f; }
done
if [ ! -s Fe.pbe-spn-kjpaw_psl.0.2.1.UPF ] || [ ! -s F.pbe-n-kjpaw_psl.0.1.UPF ]; then
  module purge; module load comp/oneapi/2023 lib/hdf5/1.14.3-oneapi-2023.0 apps/espresso/7.2-oneapi-2023.0
  [ -d pslibrary ] || git clone -q --depth 1 https://github.com/dalcorso/pslibrary.git
  awk 'NR>=482 && NR<=541' pslibrary/paw_ps_collection.job | sed 's/\$fct/pbe/g; s/\$nrel/1/g; s/\$gfun/PBE/g' > gen.sh
  bash gen.sh
  for p in Fe.pbe-spn-kjpaw_psl.0.2.1 F.pbe-n-kjpaw_psl.0.1; do [ -s $p.UPF ] || ld1.x < $p.in > $p.ld1.out; done
fi
ls -la *.UPF; grep -H "Suggested minimum cutoff" *.UPF

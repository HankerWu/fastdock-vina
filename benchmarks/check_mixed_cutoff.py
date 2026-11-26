"""Check deterministic mixed search routing and common native ranking on 1IEP."""
import argparse
from pathlib import Path
import re
import subprocess
from tempfile import TemporaryDirectory


BASIC = Path(__file__).resolve().parents[1] / 'examples/basic'


def trajectories(log):
    return re.findall(r'^Search trajectory: .*$', log, re.M)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('binary', type=Path)
    binary = str(p.parse_args().binary.resolve())
    with TemporaryDirectory() as directory:
        tmp = Path(directory)
        base = [binary, '--receptor', BASIC / '1iep_receptor.pdbqt', '--config', BASIC / '1iep_receptor.box.txt',
                '--ligand', BASIC / '1iep_ligand.pdbqt']

        def dock(name, options=(), workers=1):
            pool = tmp / (name+'.pdbqt')
            log = subprocess.check_output(base+['--seed', '511', '--cpu', str(workers), '--exhaustiveness', '3',
                    '--max_evals', '5000', '--num_modes', '6', '--energy_range', '100',
                    '--search_counters', '--out', pool, *options], text=True)
            return log, pool.read_text()

        native, _ = dock('native', ['--pose_selection', 'late'])
        smooth, _ = dock('smooth', ['--pose_selection', 'late', '--smooth_cutoff'])
        mixed, poses = dock('mixed', ['--mixed_cutoff_search'])
        mt, mt_poses = dock('threads', ['--pose_selection', 'late', '--mixed_cutoff_search'], workers=2)
        n, s = trajectories(native), trajectories(smooth)
        assert len(n) == len(s) == 3 and n[1] != s[1]
        assert trajectories(mixed) == [n[0], s[1], n[2]] == trajectories(mt)
        assert poses == mt_poses
        assert 'native=2 smooth=1' in mixed

        # With both Gaussian weights zero, no ordinary 1IEP pair term extends into
        # the 7-8 Å shell. Mixed and native Late must then agree pose for pose.
        short = ['--weight_gauss1', '0', '--weight_gauss2', '0']
        base_log, baseline = dock('short', ['--pose_selection', 'late', *short])
        other_log, other = dock('short-mixed', ['--mixed_cutoff_search', *short])
        assert baseline == other and trajectories(base_log) == trajectories(other_log)

        # Compare in the same chemical topology. PDBQT rounds coordinates to 0.001 Å.
        for i, model in enumerate(re.findall(r'^MODEL .*?^ENDMDL$', poses, re.M | re.S)):
            path = tmp / f'model{i}.pdbqt'
            path.write_text(model+'\n')
            saved = float(re.search(r'^REMARK INTER:\s+(\S+)', model, re.M)[1])
            scores = []
            for options in ([], ['--smooth_cutoff']):
                result = subprocess.check_output(base+['--pose_coordinates', path,
                    '--score_only', '--unbound_energy', '0', *options], text=True)
                scores.append(float(re.search(r'Ligand - Receptor\s+:\s+(\S+)', result)[1]))
            assert abs(scores[0]-saved) < 0.06 and abs(scores[1]-saved) > 0.5, (saved, scores)
        assert i >= 1
        # An oxygen/nitrogen-bearing ligand follows a carbon macrocycle to force new maps in batch.
        ligands = [BASIC.parents[1] / 'benchmarks/fixtures/cyclododecane.pdbqt', BASIC / '1iep_ligand.pdbqt']
        options = ['--cpu', '1', '--seed', '511', '--exhaustiveness', '3', '--max_evals', '1500',
                   '--energy_range', '100', '--num_modes', '4', '--mixed_cutoff_search']
        receptor_cmd = base[:-2]
        subprocess.run(receptor_cmd+options+['--batch', *ligands, '--dir', tmp], check=True, stdout=subprocess.PIPE)
        for ligand in ligands:
            single = tmp / (ligand.stem+'.single.pdbqt')
            subprocess.run(receptor_cmd+options+['--ligand', ligand, '--out', single], check=True, stdout=subprocess.PIPE)
            assert single.read_bytes() == (tmp / (ligand.stem+'_out.pdbqt')).read_bytes()

        for option in (['--smooth_cutoff'], ['--pose_selection', 'native'], ['--exhaustiveness', '1'],
                       ['--score_only'], ['--scoring', 'vinardo'], ['--interaction_cutoff', '1.0']):
            result = subprocess.run(base+['--cpu', '1', '--mixed_cutoff_search', '--out', tmp/'invalid.pdbqt', *option],
                                    capture_output=True, text=True)
            assert result.returncode and 'mixed_cutoff_search' in result.stderr, option
        print('Trajectory routing, native scores, zero-tail control, threads, batch and CLI constraints passed')


if __name__ == '__main__':
    main()

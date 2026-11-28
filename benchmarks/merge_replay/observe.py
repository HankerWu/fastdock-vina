"""Instrument only the post-worker merge of a clean historical checkout."""
from pathlib import Path
import shutil
import sys

here = Path(__file__).resolve().parent
lib = Path(sys.argv[1]).resolve() / "src/lib"
source = lib / "parallel_mc.cpp"
# Use only the exact conformation serialization function from the existing observer.
shutil.copyfile(here.parent / "bfgs_source_replay/trace.h", lib / "bfgs_trace.h")
old = source.read_text()
before = '''	for(const auto& trajectory : poses)
		VINA_FOR_IN(i, trajectory)
			add_to_output_container(out, trajectory[i], merge_rmsd, mc.num_saved_mins);
	out.sort();'''
after = '''	std::ofstream merge_log(std::getenv("FD_MERGE_TRACE"));
	VINA_CHECK(merge_log.good() && out.empty());
	merge_log << std::setprecision(17);
	merge_log << "P\\t" << merge_rmsd << '\\t' << mc.num_saved_mins << '\\n';
	output_container all;
	std::map<std::string, sz> identifiers;
	auto members = [&](const char* label, const output_container& pool) {
		merge_log << label;
		for (const auto& entry : pool)
			merge_log << '\\t' << identifiers.at(bfgs_trace::state(entry.c));
		merge_log << '\\n';
	};
	VINA_FOR_IN(t, poses) {
		const auto& trajectory = poses[t];
		VINA_FOR_IN(i, trajectory) {
			const auto& pose = trajectory[i];
			const auto state = bfgs_trace::state(pose.c);
			// Distinct full-precision states are required for the trace IDs below.
			VINA_CHECK(identifiers.emplace(state, all.size()).second);
			merge_log << "C\\t" << t << '\\t' << i << '\\t' << pose.e << '\\t' << state;
			for (const auto& xyz : pose.coords)
				for (sz j=0; j<3; ++j) merge_log << '\\t' << xyz[j];
			merge_log << '\\n';
			all.push_back(new output_type(pose));
			add_to_output_container(out, pose, merge_rmsd, mc.num_saved_mins);
			members("M", out);
		}
	}
	out.sort();
	members("I", out);
	const auto sorted = select_distinct_poses(all, merge_rmsd, mc.num_saved_mins);
	members("G", sorted);'''
assert old.count(before) == 1
old = old.replace('#include "coords.h"', '#include "coords.h"\n#include "bfgs_trace.h"\n#include <map>')
source.write_text(old.replace(before, after))

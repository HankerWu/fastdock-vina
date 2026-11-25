/*

   Copyright (c) 2006-2010, The Scripps Research Institute

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Author: Dr. Oleg Trott <ot14@columbia.edu>, 
           The Olson Lab, 
           The Scripps Research Institute

*/

// FastDock: reuse one model per search worker and merge by trajectory index.

#include <atomic>
#include <future>

#include "parallel_mc.h"
#include "coords.h"
#include "parallel_progress.h"
#include "trajectory_archive.h"

void parallel_mc::operator()(const model& m, output_container& out, const precalculate_byatom& p, const igrid& ig,
                             const vec& corner1, const vec& corner2, rng& generator, std::function<void(double)>* progress_callback,
                             std::vector<trajectory_archive>* archives, bool keep_support,
                             std::vector<search_trajectory_counters>* counters) const {
	VINA_CHECK(num_tasks > 0 && num_threads > 0);
	parallel_progress pp(progress_callback);
	incrementable* progress = display_progress ? &pp : nullptr;
	std::vector<rng::result_type> seeds(num_tasks);
	VINA_FOR(i, num_tasks)
		seeds[i] = static_cast<rng::result_type>(random_int(0, 1000000, generator));
	if(display_progress)
		pp.init(num_tasks * mc.global_steps);

	std::vector<output_container> poses(num_tasks);
	if (archives)
		archives->resize(num_tasks);
	if (counters)
		counters->assign(num_tasks, {});
	const vec box_center = 0.5 * (corner1 + corner2);
	std::atomic<sz> next_task{0};
	const auto run_worker = [&] {
		sz task = next_task.fetch_add(1, std::memory_order_relaxed);
		if(task >= num_tasks)
			return;
		model working(m);
		for(;;) {
			rng task_generator(seeds[task]);
			boost::optional<vec> pocket_start;
			// Even 1-based trajectories use a separate placement stream.
			if (pocket_start_radius > 0 && task % 2 == 1) {
				rng placement_generator(seeds[task] ^ 0xa511e9b3u);
				pocket_start = box_center + pocket_start_radius * random_inside_sphere(placement_generator);
			}
			mc(working, poses[task], p, ig, corner1, corner2, progress, task_generator,
			   archives ? &(*archives)[task] : nullptr, counters ? &(*counters)[task] : nullptr,
			   pocket_start);
			task = next_task.fetch_add(1, std::memory_order_relaxed);
			if(task >= num_tasks)
				break;
			// The first mutation can inspect the starting geometry before evaluation updates it.
			working.reset_search_state(m);
		}
	};

	const sz workers = std::min(num_threads, num_tasks);
	std::vector<std::future<void>> background;
	background.reserve(workers - 1);
	for(sz i = 1; i < workers; ++i)
		background.push_back(std::async(std::launch::async, run_worker));
	run_worker();
	for(auto& worker : background)
		worker.get();

	if(keep_trajectory_minima) {
		sz total = 0;
		for(const auto& trajectory : poses)
			total += trajectory.size();
		out.reserve(out.size() + total);
		if (keep_support) {
			VINA_FOR(i, num_tasks) {
				auto& archive = (*archives)[i];
				archive.samples.reserve(poses[i].size() + archive.accepted.size());
				for (auto& pose : poses[i])
					archive.samples.push_back(pose.coords);
				for (auto& pose : archive.accepted)
					archive.samples.push_back(pose.coords);
			}
		}
		for(auto& trajectory : poses)
			out.transfer(out.end(), trajectory);
		return;
	}

	// Preserve the historical trajectory order and ties.
	const fl merge_rmsd = use_min_rmsd_for_merge ? mc.min_rmsd : 2;
	for(const auto& trajectory : poses)
		VINA_FOR_IN(i, trajectory)
			add_to_output_container(out, trajectory[i], merge_rmsd, mc.num_saved_mins);
	out.sort();
}

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

void parallel_mc::operator()(const model& m, output_container& out, const precalculate_byatom& p, const igrid& ig, const vec& corner1, const vec& corner2, rng& generator, std::function<void(double)>* progress_callback) const {
	VINA_CHECK(num_tasks > 0 && num_threads > 0);
	parallel_progress pp(progress_callback);
	incrementable* progress = display_progress ? &pp : nullptr;
	std::vector<rng::result_type> seeds(num_tasks);
	VINA_FOR(i, num_tasks)
		seeds[i] = static_cast<rng::result_type>(random_int(0, 1000000, generator));
	if(display_progress)
		pp.init(num_tasks * mc.global_steps);

	std::vector<output_container> poses(num_tasks);
	std::atomic<sz> next_task{0};
	const auto run_worker = [&] {
		sz task = next_task.fetch_add(1, std::memory_order_relaxed);
		if(task >= num_tasks)
			return;
		model working(m);
		for(;;) {
			rng task_generator(seeds[task]);
			mc(working, poses[task], p, ig, corner1, corner2, progress, task_generator);
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

	// Keep the historical merge distance and trajectory order, including ties.
	for(const auto& trajectory : poses)
		VINA_FOR_IN(i, trajectory)
			add_to_output_container(out, trajectory[i], 2, mc.num_saved_mins);
	out.sort();
}

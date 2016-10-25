def create_file_targets(variables, config):

    merge_config = config["merge"]

    merge_mapper = {}
    original_lane_mapper = {}
    final_samples = []
    run_lanes = []
    run_merged_lane = {}
    unmerged_samples = []
    final_samples = []
    all_samples = []

    # Create mappers and new run_lane combinations based on config
    for run in merge_config:
        run_merged_lane[str(run)] = []
        for merge in merge_config[run]:
            merge = [str(l) for l in merge]
            run_merged_lane[str(run)] += merge
            new_lane = "{run}_{lane1}-{lane2}".format(
                run=run, lane1=merge[0], lane2=merge[1])
            run_lanes.append(new_lane)
            run_lane1 = "{run}_{lane}".format(run=run, lane=merge[0])
            run_lane2 = "{run}_{lane}".format(run=run, lane=merge[1])
            original_lane_mapper[run_lane1] = [merge[0], merge[1]]
            original_lane_mapper[run_lane2] = [merge[0], merge[1]]

    # Run through all the runs, lanes and tags found in the cram folder
    for i in range(len(variables.run)):
        run, lane, tag_index = (
            variables.run[i], variables.lane[i], variables.tag_index[i])
        run_lane = str(run) + "_" + str(lane)
        original_sample_name = "{run_lane}#{tag_index}".format(
            run_lane=run_lane, tag_index=tag_index)
        all_samples.append(original_sample_name)

        if lane not in run_merged_lane[run]:
            if run_lane not in run_lanes:
                run_lanes.append(run_lane)
            # merge_mapper[sample_name] = [original_sample_name, None]
            unmerged_samples.append(original_sample_name)
        else:
            lane1 = original_lane_mapper[run_lane][0]
            lane2 = original_lane_mapper[run_lane][1]
            # Get the new sample name
            merged_lane = "{run}_{lane1}-{lane2}".format(
                run=run, lane1=lane1, lane2=lane2)
            sample_name = "{run_lane}#{tag_index}".format(
                run_lane=merged_lane, tag_index=tag_index)

            # Map the old sample names for the merge
            original_sample1 = config["pattern"].format(
                run=run, lane=lane1, tag_index=tag_index)
            original_sample2 = config["pattern"].format(
                run=run, lane=lane2, tag_index=tag_index)
            merge_mapper[sample_name] = [original_sample1, original_sample2]

    # Look for missing files
    final_merge_mapper = {}
    for sample in merge_mapper:
        if (merge_mapper[sample][0] in all_samples and
              merge_mapper[sample][1] in all_samples):
            final_merge_mapper[sample] = merge_mapper[sample]
        else:
            print("Could not find both samples for merge: %s" %
                  ", ".join(merge_mapper[sample]))

    final_samples = list(final_merge_mapper.keys()) + unmerged_samples

    return final_samples, final_merge_mapper



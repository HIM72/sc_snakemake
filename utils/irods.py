from __future__ import print_function

import six

import hashlib
import os
import sys
import copy
import json
import subprocess
from pprint import pprint
import pandas as pd


class Collection(object):

    def __init__(self, study_id):
        self.study_id = study_id

    @property
    def base_query(self):
        return {"collection": "/seq",
                "avus": [{"attribute": "target", "value": "1"},
                         {"attribute": "study_id", "value": "%s" % self.study_id}]}
                        #{"attribute": "type", "value": "cram"},

    def _convert_avus(self, avus):
        return {a['attribute']: a['value'] for a in avus}

    def create_query(self, baton_command, query=None, **avus):
        if not query:
            query = copy.deepcopy(self.base_query)
        for key, value in six.iteritems(avus):
            if value is not None:
                query["avus"].append({"attribute": key, "value": str(value)})
        return "jq -n '%s' | %s" % (json.dumps(query), baton_command)

    def execute_baton(self, baton_command, query=None, **avus):
        command = self.create_query(baton_command, query=query, **avus)
        print("Executing baton command: %s" % baton_command)
        try:
            result = subprocess.getoutput(command)
            try:
                return json.loads(result)
            except (TypeError, ValueError):
                return
        except subprocess.CalledProcessError:
            print("Could not execute: '%s', '%s', '%s'" %
                  (baton_command, str(query), str(avus)))

    def move_from_irods(self, out_folder=None, id_run=None, lane=None,
                        tag_index=None):
        if not out_folder:
            out_folder = '/nfs/team205/.scapi/%s/cram' % str(self.study_id)
        if not os.path.isdir(out_folder):
            os.makedirs(out_folder)
        data_objects = self.execute_baton(
            'baton-metaquery --avu', id_run=id_run, lane=lane,
            tag_index=tag_index)
        get_query = copy.deepcopy(data_objects)
        # get_query = self.iget_targets(out_folder, id_run=id_run, lane=lane)
        for index, item in list(enumerate(get_query)):
            avus = self._convert_avus(item.pop('avus'))
            filename = item['data_object']
            if 'phix' not in filename:
                item["directory"] = out_folder
                full_path = os.path.join(out_folder, filename)
                data_objects[index]['location'] = full_path
                if not os.path.isfile(full_path):
                    result = self.execute_baton('baton-get --save', query=item)
                else:
                    if avus["md5"]:
                        checksum = hashlib.md5(open(full_path, 'rb').read()) \
                            .hexdigest()
                        if checksum != avus["md5"]:
                            os.remove(full_path)
                            result = self.execute_baton('baton-get --save',
                                                        query=item)

        return data_objects

    def iget_targets(self, out_folder=None, id_run=None, lane=None):
        if not out_folder:
            out_folder = '/nfs/team205/.scapi/%s/cram' % str(self.study_id)
        if not os.path.isdir(out_folder):
            os.makedirs(out_folder)
        # query = copy.deepcopy(self.base_query)
        # if id_run:
        #     query["avus"].append({"attribute": "id_run", "value": str(id_run)})
        # if lane:
        #     query["avus"].append({"attribute": "lane", "value": str(lane)})
        data_objects = self.execute_baton('baton-metaquery',
                                          id_run=id_run, lane=lane)
        print(data_objects)
        get_query = copy.deepcopy(data_objects)
        for item in get_query:
            item['directory'] = out_folder

        return get_query

    def check_metadata(self):
        # Get the metadata of the first object
        data_objects = self.execute_baton('baton-metaquery')
        metadata = self.execute_baton('baton-list --avu', query=data_objects[0])

        if metadata:
            return self._convert_avus(metadata['avus'])
        else:
            return

    def get_collection_metadata(self, id_run=None, lane=None, tag_index=None):
        data_objects = self.execute_baton(
            'baton-metaquery --avu', id_run=id_run,
            lane=lane, tag_index=tag_index)
        if data_objects:
            for data_object in data_objects:
                data_object['avus'] = self._convert_avus(data_object['avus'])

            return data_objects
        else:
            return []


if __name__ == '__main__':
    c = Collection(study_id=18482)
    pprint(c.get_collection_metadata()[1])

/*
  Copyright (c) 2003-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "gtt.h"
#include "core/array.h"
#include "core/array2dim.h"
#include "core/bitpackarray.h"
#include "core/bitpackstring.h"
#include "core/bittab.h"
#include "core/bsearch.h"
#include "core/countingsort.h"
#include "core/cstr_table.h"
#include "core/disc_distri.h"
#include "core/dlist.h"
#include "core/dynbittab.h"
#include "core/getbasename.h"
#include "core/grep.h"
#include "core/hashmap.h"
#include "core/hashtable.h"
#include "core/interval_tree.h"
#include "core/queue.h"
#include "core/splitter.h"
#include "core/tokenizer.h"
#include "extended/alignment.h"
#include "extended/evaluator.h"
#include "extended/genome_feature.h"
#include "extended/genome_node_iterator.h"
#include "extended/gff3_escaping.h"
#include "extended/hmm.h"
#include "extended/luaserialize.h"
#include "extended/splicedseq.h"
#include "extended/string_matching.h"
#include "extended/tag_value_map.h"
#include "extended/union_find.h"
#include "extended/redblack.h"
#include "tools/gt_bioseq.h"
#include "tools/gt_cds.h"
#include "tools/gt_chseqids.h"
#include "tools/gt_clean.h"
#include "tools/gt_csa.h"
#include "tools/gt_dev.h"
#include "tools/gt_eval.h"
#include "tools/gt_exercise.h"
#include "tools/gt_extractfeat.h"
#include "tools/gt_extractseq.h"
#include "tools/gt_filter.h"
#include "tools/gt_fingerprint.h"
#include "tools/gt_gff3.h"
#include "tools/gt_gff3validator.h"
#include "tools/gt_gff3_to_gtf.h"
#include "tools/gt_gtf_to_gff3.h"
#include "tools/gt_ltrharvest.h"
#include "tools/gt_matchingstatistics.h"
#include "tools/gt_merge.h"
#include "tools/gt_mgth.h"
#include "tools/gt_mkfmindex.h"
#include "tools/gt_mmapandread.h"
#include "tools/gt_mutate.h"
#include "tools/gt_packedindex.h"
#include "tools/gt_prebwt.h"
#include "tools/gt_seqfilter.h"
#include "tools/gt_sequniq.h"
#include "tools/gt_shredder.h"
#include "tools/gt_splitfasta.h"
#include "tools/gt_splicesiteinfo.h"
#include "tools/gt_stat.h"
#include "tools/gt_suffixerator.h"
#include "tools/gt_tagerator.h"
#include "tools/gt_template.h"
#include "tools/gt_uniq.h"
#include "tools/gt_uniquesub.h"

#ifndef WITHOUT_CAIRO
#include "annotationsketch/block.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/feature_index.h"
#include "annotationsketch/gt_sketch.h"
#include "annotationsketch/image_info.h"
#include "annotationsketch/track.h"
#include "annotationsketch/recmap.h"
#include "annotationsketch/style.h"
#endif

Toolbox* gtt_tools(void)
{
  Toolbox *tools = toolbox_new();

  /* add tools */
  toolbox_add_tool(tools, "bioseq", gt_bioseq());
  toolbox_add_tool(tools, "cds", gt_cds());
  toolbox_add(tools, "chseqids", gt_chseqids);
  toolbox_add(tools, "clean", gt_clean);
  toolbox_add_tool(tools, "csa", gt_csa());
  toolbox_add_tool(tools, "dev", gt_dev());
  toolbox_add(tools, "eval", gt_eval);
  toolbox_add_tool(tools, "exercise", gt_exercise());
  toolbox_add_tool(tools, "extractfeat", gt_extractfeat());
  toolbox_add_tool(tools, "extractseq", gt_extractseq());
  toolbox_add_tool(tools, "filter", gt_filter());
  toolbox_add_tool(tools, "fingerprint", gt_fingerprint());
  toolbox_add_tool(tools, "gff3", gt_gff3());
  toolbox_add_tool(tools, "gff3validator", gt_gff3validator());
  toolbox_add(tools, "gff3_to_gtf", gt_gff3_to_gtf);
  toolbox_add(tools, "gtf_to_gff3", gt_gtf_to_gff3);
  toolbox_add(tools, "ltrharvest", gt_ltrharvest);
  toolbox_add(tools, "matstat", gt_matchingstatistics);
  toolbox_add(tools, "merge", gt_merge);
  toolbox_add(tools, "mgth", gt_mgth);
  toolbox_add(tools, "mmapandread", gt_mmapandread);
  toolbox_add_tool(tools, "mutate", gt_mutate());
  toolbox_add(tools, "mkfmindex", gt_mkfmindex);
  toolbox_add_tool(tools, "packedindex", gt_packedindex());
  toolbox_add_tool(tools, "prebwt", gt_prebwt());
  toolbox_add_tool(tools, "seqfilter", gt_seqfilter());
  toolbox_add_tool(tools, "sequniq", gt_sequniq());
  toolbox_add_tool(tools, "shredder", gt_shredder());
  toolbox_add_tool(tools, "splitfasta", gt_splitfasta());
  toolbox_add(tools, "splicesiteinfo", gt_splicesiteinfo);
  toolbox_add(tools, "stat", gt_stat);
  toolbox_add(tools, "suffixerator", gt_suffixerator);
  toolbox_add_tool(tools, "tagerator", gt_tagerator());
  toolbox_add_tool(tools, "template", gt_template());
  toolbox_add(tools, "uniq", gt_uniq);
  toolbox_add(tools, "uniquesub", gt_uniquesub);
#ifndef WITHOUT_CAIRO
  toolbox_add(tools, "sketch", gt_sketch);
#endif

  return tools;
}

Hashmap* gtt_unit_tests(void)
{
  Hashmap *unit_tests = hashmap_new(HASH_STRING, NULL, NULL);

  /* add unit tests */
  hashmap_add(unit_tests, "alignment class", alignment_unit_test);
  hashmap_add(unit_tests, "array class", gt_array_unit_test);
  hashmap_add(unit_tests, "array example", gt_array_example);
  hashmap_add(unit_tests, "array2dim example", array2dim_example);
  hashmap_add(unit_tests, "bit pack array class", gt_bitpackarray_unit_test);
  hashmap_add(unit_tests, "bit pack string module", bitPackString_unit_test);
  hashmap_add(unit_tests, "bittab class", gt_bittab_unit_test);
  hashmap_add(unit_tests, "bittab example", gt_bittab_example);
  hashmap_add(unit_tests, "bsearch module", bsearch_unit_test);
  hashmap_add(unit_tests, "countingsort module", countingsort_unit_test);
  hashmap_add(unit_tests, "cstr table class", cstr_table_unit_test);
  hashmap_add(unit_tests, "disc distri class", disc_distri_unit_test);
  hashmap_add(unit_tests, "dlist class", gt_dlist_unit_test);
  hashmap_add(unit_tests, "dlist example", gt_dlist_example);
  hashmap_add(unit_tests, "dynamic bittab class", dynbittab_unit_test);
  hashmap_add(unit_tests, "evaluator class", evaluator_unit_test);
  hashmap_add(unit_tests, "genome feature class", gt_genome_feature_unit_test);
  hashmap_add(unit_tests, "genome node iterator example",
              gt_genome_node_iterator_example);
  hashmap_add(unit_tests, "getbasename module", getbasename_unit_test);
  hashmap_add(unit_tests, "gff3 escaping module", gff3_escaping_unit_test);
  hashmap_add(unit_tests, "grep module", grep_unit_test);
  hashmap_add(unit_tests, "hashmap class", hashmap_unit_test);
  hashmap_add(unit_tests, "hashtable class", hashtable_unit_test);
  hashmap_add(unit_tests, "hmm class", hmm_unit_test);
  hashmap_add(unit_tests, "interval tree class", interval_tree_unit_test);
  hashmap_add(unit_tests, "Lua serializer module", lua_serializer_unit_test);
  hashmap_add(unit_tests, "queue class", queue_unit_test);
  hashmap_add(unit_tests, "range class", gt_range_unit_test);
  hashmap_add(unit_tests, "red-black tree class", rbt_unit_test);
  hashmap_add(unit_tests, "safearith module", safearith_unit_test);
  hashmap_add(unit_tests, "safearith example", safearith_example);
  hashmap_add(unit_tests, "splicedseq class", splicedseq_unit_test);
  hashmap_add(unit_tests, "splitter class", splitter_unit_test);
  hashmap_add(unit_tests, "string class", gt_str_unit_test);
  hashmap_add(unit_tests, "string matching module", string_matching_unit_test);
  hashmap_add(unit_tests, "tag value map example", tag_value_map_example);
  hashmap_add(unit_tests, "tokenizer class", tokenizer_unit_test);
  hashmap_add(unit_tests, "union find class", union_find_unit_test);
#ifndef WITHOUT_CAIRO
  hashmap_add(unit_tests, "block class", gt_block_unit_test);
  hashmap_add(unit_tests, "style class", gt_style_unit_test);
  hashmap_add(unit_tests, "diagram class", gt_diagram_unit_test);
  hashmap_add(unit_tests, "element class", gt_element_unit_test);
  hashmap_add(unit_tests, "feature index class", gt_feature_index_unit_test);
  hashmap_add(unit_tests, "imageinfo class", gt_image_info_unit_test);
  hashmap_add(unit_tests, "line class", gt_line_unit_test);
  hashmap_add(unit_tests, "track class", gt_track_unit_test);
#endif

  return unit_tests;
}

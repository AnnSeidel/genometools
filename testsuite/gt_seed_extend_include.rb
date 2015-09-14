def build_encseq(indexname, sequencefile)
  return "#{$bin}gt encseq encode -des no -sds no -md5 no " +
    "-indexname " + indexname + " " + sequencefile
end

Name "gt seed_extend mirror, check k-mers and seed pairs"
Keywords "gt_seed_extend seedpair kmer polysequence xdrop extend"
Test do
  run_test build_encseq("small_poly", "#{$testdata}small_poly.fas")
  run_test "#{$bin}gt seed_extend -seedlength 10 -mirror -debug-kmer " +
           "-debug-seedpair small_poly"
  run "cmp -s #{last_stdout} #{$testdata}seedextend1.out"
  run_test "#{$bin}gt seed_extend -seedlength 10 -mirror -extendxdrop 97 " +
           "-l 10 -mincoverage 11 small_poly"
  run "cmp -s #{last_stdout} #{$testdata}seedextend3.out"
end

Name "gt seed_extend memlimit, use wildcard containing reads"
Keywords "gt_seed_extend seedpair at1MB memlimit maxfreq verbose"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test "#{$bin}gt seed_extend -verify -debug-seedpair -memlimit 10MB at1MB"
  run "cmp -s #{last_stdout} #{$testdata}seedextend2.out"
  run_test "#{$bin}gt seed_extend -v -maxfreq 5 at1MB"
  grep last_stdout, /...found and sorted 582230 k-mers/
  grep last_stdout, /...collected and sorted 68577 seed pairs/
end

Name "gt seed_extend failure"
Keywords "gt_seed_extend fail"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test build_encseq("foo", "#{$testdata}foo.fas")
  run_test "#{$bin}gt seed_extend -seedlength 15 at1MB", :retval => 1
  grep last_stderr, /integer <= 14 if the sequences contain wildcards/
  run_test "#{$bin}gt seed_extend -seedlength 10 foo", :retval => 1
  grep last_stderr, /integer <= 8 \(length of longest sequence\)/
  run_test "#{$bin}gt seed_extend -maxfreq 1 at1MB", :retval => 1
  grep last_stderr, /option "-maxfreq" must be >= 2 to find matching k-mers/
  run_test "#{$bin}gt seed_extend -t 2 at1MB", :retval => 1
  grep last_stderr, /option "-t" must be >= 3 to find matching k-mers/
  run_test "#{$bin}gt seed_extend -memlimit 0MB at1MB", :retval => 1
  grep last_stderr, /argument to option "-memlimit" must be at least 1MB/
  run_test "#{$bin}gt seed_extend -memlimit 1MB at1MB", :retval => 1
  grep last_stderr, /option -memlimit too strict: need at least 10MB/
  run_test "#{$bin}gt seed_extend -memlimit 1KB at1MB", :retval => 1
  grep last_stderr, /integer argument followed by one of the keywords MB and GB/
  run_test "#{$bin}gt seed_extend -extendgreedy -history 65 -benchmark at1MB",
    :retval => 1
  grep last_stderr, /argument to option "-history" must be an integer <= 64/
  run_test "#{$bin}gt seed_extend -percmathistory 140 -extendgreedy -v at1MB",
    :retval => 1
  grep last_stderr, /option "-percmathistory" must be an integer <= 100/
  run_test "#{$bin}gt seed_extend -extendgreedy -cam invalidlongcamstring " +
    "at1MB", :retval => 1
  grep last_stderr, /illegal parameter for option -cam/
  run_test "#{$bin}gt seed_extend -v at1MB at1MB at1MB", :retval => 1
  grep last_stderr, /too many arguments/
  run_test "#{$bin}gt seed_extend -benchmark", :retval => 1
  grep last_stderr, /at least one encseq index name must be specified/
  run_test "#{$bin}gt seed_extend -extendgreedy -seedlength 7 at1MB foo",
    :retval => 1
  grep last_stderr, /comparison of two encseqs not implemented/
end

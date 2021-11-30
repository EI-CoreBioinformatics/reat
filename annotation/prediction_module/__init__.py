def add_classification_parser_parameters(parser):
    parser.add_argument('--min_pct_cds_fraction', type=float, default=0.5,
                        help='Transcript requirement for minimum faction of transcript covered by CDS')
    parser.add_argument('--max_tp_utr_complete', type=int, default=1,
                        help='Transcript requirement for maximum complete 3\' UTRs')
    parser.add_argument('--max_tp_utr', type=int, default=2,
                        help='Transcript requirement for maximum number of 3\' UTRs')
    parser.add_argument('--min_tp_utr', type=int, default=1,
                        help='Transcript requirement for minumum number of 3\' UTRs')
    parser.add_argument('--max_fp_utr_complete', type=int, default=2,
                        help='Transcript requirement for maximum number of complete 5\' UTRs')
    parser.add_argument('--max_fp_utr', type=int, default=3,
                        help='Transcript requirement for maximum number of 5\' UTRs')
    parser.add_argument('--min_fp_utr', type=int, default=1,
                        help='Transcript requirement for minumum number of 5\' UTRs')
    parser.add_argument('--query_start_hard_filter_distance', type=int, default=10,
                        help='If query hit starts after this value, the transcript cannot belong to the Gold category')
    parser.add_argument('--query_start_score', type=int, default=5,
                        help='Maximum score for query start distance')
    parser.add_argument('--query_start_scoring_distance', type=int, default=30,
                        help='Hits with query start distance lower than this parameter start receiving scoring points')
    parser.add_argument('--query_end_hard_filter_distance', type=int, default=10,
                        help='If query hit ends after this value, the transcript cannot belong to the Gold category')
    parser.add_argument('--query_end_score', type=int, default=5,
                        help='Maximum score for query end distance')
    parser.add_argument('--query_end_scoring_distance', type=int, default=30,
                        help='Hits with query end distance lower than this parameter start receiving scoring points')
    parser.add_argument('--target_start_hard_filter_distance', type=int, default=10,
                        help='If target hit starts after this value, the transcript cannot belong to the Gold category')
    parser.add_argument('--target_start_score', type=int, default=5,
                        help='Maximum score for target start distance')
    parser.add_argument('--target_start_scoring_distance', type=int, default=30,
                        help='Hits with target start distance lower than this parameter start receiving scoring points')
    parser.add_argument('--target_end_hard_filter_distance', type=int, default=10,
                        help='If target hit ends after this value, the transcript cannot belong to the Gold category')
    parser.add_argument('--target_end_score', type=int, default=5,
                        help='Maximum score for target end distance')
    parser.add_argument('--target_end_scoring_distance', type=int, default=30,
                        help='Hits with target end distance lower than this parameter start receiving scoring points')
    parser.add_argument('--min_query_coverage_hard_filter', type=int, default=90,
                        help='Minimum percentage of query covered to classify a hit as Gold')
    parser.add_argument('--min_query_coverage_score', type=int, default=5,
                        help='Maximum score for query percentage coverage')
    parser.add_argument('--min_query_coverage_scoring_percentage', type=int, default=30,
                        help='Queries covered over this percentage value start receiving scoring points')
    parser.add_argument('--min_target_coverage_hard_filter', type=int, default=90,
                        help='Minimum percentage of target covered to classify a hit as Gold')
    parser.add_argument('--min_target_coverage_score', type=int, default=5,
                        help='Maximum score for target percentage coverage')
    parser.add_argument('--min_target_coverage_scoring_percentage', type=int, default=30,
                        help='Targets covered over this percentage value start receiving scoring points')
    parser.add_argument('--max_single_gap_hard_filter', type=int, default=20,
                        help='Any hits containing gaps larger than this parameter cannot be classified as Gold')
    parser.add_argument('--max_single_gap_score', type=int, default=5,
                        help='Maximum score for hits with single gaps smaller than --max_single_gap_scoring_length')
    parser.add_argument('--max_single_gap_scoring_length', type=int, default=30,
                        help='Hits containing gaps smaller than this parameter start receiving scoring points')
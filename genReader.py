def generator_to_dicts(S):
    new_sp = {si: {"list": lists[1]} for si, lists in S.sp.items()}
    for si, si_info in new_sp.items():
        si_info["list_rank"] = dict()
        for tie_idx, tie in enumerate(si_info["list"]):
            for pj_idx, pj in enumerate(tie):
                si_info["list_rank"][pj] = (tie_idx, pj_idx)

        si_info["list_len"] = len(si_info["list"])
        si_info["head_idx"] = 0

    new_sp_no_tie_deletions = {si: set(lists[0]) for si, lists in S.sp.items()}

    new_lp = {
        lk: {
            "cap": lk_details[0],
            "list": lk_details[-1],
            "list_len": len(lk_details[-1]),
            "tail_idx": len(lk_details[-1]) - 1,
            "projects": set(lk_details[1]),
        }
        for lk, lk_details in S.lp.items()
    }

    new_plc = {
        pj: {
            "cap": pj_details[0],
            "lec": pj_details[1],
        }
        for pj, pj_details in S.plc.items()
    }
    for pj, pj_details in new_plc.items():
        lk = new_plc[pj]["lec"]
        accessible_students = S.plc[pj][2]

        pj_details["list"] = []
        for tie in new_lp[lk]["list"].copy():
            new_tie = []
            for si in tie:
                if si in accessible_students:
                    new_tie.append(si)
            if new_tie:
                pj_details["list"].append(new_tie)

        pj_details["list_len"] = len(pj_details["list"])
        pj_details["tail_idx"] = pj_details["list_len"] - 1

    return new_sp, new_sp_no_tie_deletions, new_plc, new_lp


def generator_to_bruteforcer_dicts(S):
    dict_sp, new_sp_no_tie_deletions, dict_plc, dict_lp = generator_to_dicts(S)

    new_sp = {
        si: [info["list_len"], info["list"], info["list_rank"], info["head_idx"]]
        for si, info in dict_sp.items()
    }

    new_plc = {
        pj: [info["lec"], info["cap"], False, [], info["list_len"], info["tail_idx"]]
        for pj, info in dict_plc.items()
    }

    new_lp = {
        lk: [
            info["cap"],
            info["list"],
            False,
            dict(),
            info["list_len"],
            info["tail_idx"],
        ]
        for lk, info in dict_lp.items()
    }

    new_sp_no_tie = {si: lists[0].copy() for si, lists in S.sp.items()}

    new_lp_rank = {lk: dict() for lk in S.lp}
    for lk, lk_ranking in new_lp_rank.items():
        lk_list = S.lp[lk][-1]
        for tie_idx, tie in enumerate(lk_list):
            for si in tie:
                lk_ranking[si] = tie_idx + 1

    new_proj_rank = {pj: dict() for pj in S.plc}
    for pj, pj_ranking in new_proj_rank.items():
        lk = S.plc[pj][1]
        accessible_students = S.plc[pj][2]

        for si, si_rank in new_lp_rank[lk].items():
            if si in accessible_students:
                pj_ranking[si] = si_rank

        if pj_ranking.values():
            min_rank = min(pj_ranking.values())
            for si, si_rank in pj_ranking.items():
                pj_ranking[si] = si_rank + 1 - min_rank

    return (
        new_sp,
        new_sp_no_tie,
        new_sp_no_tie_deletions,
        new_plc,
        new_lp,
        new_lp_rank,
        new_proj_rank,
    )

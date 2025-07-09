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

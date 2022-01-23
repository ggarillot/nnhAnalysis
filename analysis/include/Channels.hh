#pragma once

#include <map>

#include <TColor.h>

enum Channel
{
    SIGNAL,
    OTHER_H_WW,
    OTHER_H_NOT_WW,
    OTHER_H_BB,
    OTHER_H_NOT_BB,
    F2_L,
    F2_H,
    F4_L,
    F4_H,
    F4_SL,
    BKG
};

struct ChannelInfo
{
    std::string displayName = "";
    Color_t     color = 0;
};

const std::map<Channel, ChannelInfo> CHANNEL_PARAMS = {{SIGNAL, {"Signal", 1}},
                                                       {OTHER_H_WW, {"Other higgs->WW*", 401}},
                                                       {OTHER_H_NOT_WW, {"Other higgs", 417}},
                                                       {OTHER_H_BB, {"Other higgs->bb", 401}},
                                                       {OTHER_H_NOT_BB, {"Other higgs", 417}},
                                                       {F2_L, {"2 fermions leptonic", 435}},
                                                       {F2_H, {"2 fermions hadronic", 633}},
                                                       {F4_L, {"4 fermions leptonic", 433}},
                                                       {F4_H, {"4 fermions hadronic", 625}},
                                                       {F4_SL, {"4 fermions semileptonic", 591}},
                                                       {BKG, {"Background", 63}}};

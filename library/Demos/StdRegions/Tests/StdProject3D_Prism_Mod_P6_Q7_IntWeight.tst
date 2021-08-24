<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject3D Prism Modified basis P=6 Q=7 with alternative weights</description>
    <executable>StdProject</executable>
    <parameters>-s prism -b Modified_A Modified_A Modified_B -o 6 6 6 -p 7 7 7 -P GaussLobattoLegendre GaussLobattoLegendre GaussLobattoLegendre</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-11">2.08237e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">4.65228e-12</value>
        </metric>
        <metric type="Regex" id="3">
            <regex>^Integral error: ([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
            <matches>
                <match>
                    <field id="0" tolerance="1e-11">7.38076e-13</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>



<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>StdProject3D Pyramid Modified basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s pyramid -b Modified_A Modified_A ModifiedPyr_C -o 6 6 6 -p 7 7 7 -P GaussLobattoLegendre GaussLobattoLegendre GaussLobattoLegendre</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-11">8.8186e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">2.41681e-12</value>
        </metric>
        <metric type="Regex" id="3">
            <regex>^Integral error: ([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
            <matches>
                <match>
                    <field id="0" tolerance="1e-11">2.49289e-12</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>

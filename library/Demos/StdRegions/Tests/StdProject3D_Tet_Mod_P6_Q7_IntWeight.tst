<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>StdProject3D Tetrahedron Modified basis P=6 Q=7, alternative quadrature weights</description>
    <executable>StdProject</executable>
    <parameters>-s tetrahedron -b Modified_A Modified_B Modified_C -o 6 6 6 -p 7 7 7 -P GaussLobattoLegendre GaussLobattoLegendre GaussLobattoLegendre</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">2.89257e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="2e-12">7.77156e-13</value>
        </metric>
        <metric type="Regex" id="3">
            <regex>^Integral error: ([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
            <matches>
                <match>
                    <field id="0" tolerance="1e-11">2.15439e-13</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>

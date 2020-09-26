<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>StdProject2D Triangle Modified basis P=6 Q=7, alternative quadrature weights</description>
    <executable>StdProject</executable>
    <parameters>-s triangle -b Modified_A Modified_B -o 6 6 -p 7 7 -P GaussLobattoLegendre GaussLobattoLegendre</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">3.24325e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">3.86358e-14</value>
        </metric>
        <metric type="Regex" id="3">
            <regex>^Integral error: ([+-]?\d.+\d|-?\d|[+-]?nan|[+-]?inf).*</regex>
            <matches>
                <match>
                    <field id="0" tolerance="1e-11">3.99614e-12</field>
                </match>
            </matches>
        </metric>
    </metrics>
</test>

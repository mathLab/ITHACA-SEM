<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LocProject_Diff2D Tri Mod. Basis P=6, Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s triangle -b Modified_A Modified_B -o 6 6 -p 7 7 -c 0.0 0.0 1.0 1.0 0.5 1.0 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">7.02226e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="2e-12">1.13687e-12</value>
        </metric>
    </metrics>
</test>



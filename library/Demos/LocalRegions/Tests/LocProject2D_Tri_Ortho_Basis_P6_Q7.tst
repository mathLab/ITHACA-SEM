<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project2D Triangle Orthogonal basis P=6 Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s triangle -b Ortho_A Ortho_B -o 6 6 -p 7 7 -c 0.0 0.0 1.0 1.0 0.5 1.0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">1.25467e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">6.21725e-14</value>
        </metric>
    </metrics>
</test>



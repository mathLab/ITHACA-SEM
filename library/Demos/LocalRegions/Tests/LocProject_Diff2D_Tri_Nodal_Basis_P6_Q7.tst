<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LocProject_Diff2D Tri Nodal Basis P=6, Q=7</description>
    <executable>LocProject_Diff2D</executable>
    <parameters>3 11 11 6 6 7 7  0.0 0.0 1.0 1.0 0.5 1.0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">2.99203e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">6.11067e-13</value>
        </metric>
    </metrics>
</test>



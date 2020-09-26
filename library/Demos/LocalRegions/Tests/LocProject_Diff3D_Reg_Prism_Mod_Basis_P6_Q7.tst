<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LocProject_Diff3D Reg. Prism Modified Basis, P=6, Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s prism -b Modified_A Modified_A Modified_B -o 6 6 6 -p 7 7 7 -c 0.0 0.0 0.0 1.0 0.0 0.0 1.0 1.0 0.0 0.0 1.0 0.0 0.5 0.0 1.0 0.5 1.0 1.0 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-9">2.41945e-11</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-9">1.03926e-10</value>
        </metric>
    </metrics>
</test>



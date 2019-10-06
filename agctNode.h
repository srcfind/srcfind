//
// Created by ttbond on 18-10-17.
//

#ifndef SRCFIND_AGCTNODE_H
#define SRCFIND_AGCTNODE_H


class agctNode {
public:
    agctNode();
    inline agctNode *enterSon(int id);
    agctNode *son[4];
    long long num;
};

inline agctNode *agctNode::enterSon(int id){
    if(id<0||id>3){
        printf("wrong Id in agctNode::enterSon(int id) value:%d\n",id);
        exit(-1);
    }
    if(son[id]==NULL){
        son[id]=new agctNode();
    }
    return son[id];
}
#endif //SRCFIND_AGCTNODE_H
